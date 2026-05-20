#!/usr/bin/env python3
"""
Space-filling street network model.

Model 1 (basic): pick a segment with probability proportional to length, split
  at midpoint, connect midpoint to the closest visible intersection.

Model 2 (biased): same as Model 1 but only splits segments whose half-length²
  exceeds area_coeff·exp(-area_scale/r) (larger blocks near the periphery), and
  rejects connections that would create an intersection angle below min_angle.

Both start from a unit square divided by 4 lines (horizontal, vertical, two
diagonals) crossing at the centre.

Performance: the hot inner loops (visibility, angle check) are JIT-compiled
with numba when available; otherwise pure-numpy fallbacks are used. Use
`run_ensemble` to run many seeds in parallel.

Requires: numpy, scipy, shapely, geopandas. Optional: numba.
"""

from __future__ import annotations

import os
from collections import defaultdict
from multiprocessing import Pool

import geopandas as gpd
import numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import LineString

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return lambda f: f

MIN_LENGTH = 0.05
CENTER_X, CENTER_Y = 0.5, 0.5
EPS = 1e-9


# --------------------------------------------------------------------------- #
# Hot inner loops (JIT-compiled when numba is available)                      #
# --------------------------------------------------------------------------- #

@njit(cache=True)
def _segment_blocks(px, py, qx, qy, seg_arr):
    """True if any segment in seg_arr properly crosses p→q.

    Segments sharing an endpoint with p or q are treated as non-blocking.
    """
    for i in range(seg_arr.shape[0]):
        ax = seg_arr[i, 0, 0]; ay = seg_arr[i, 0, 1]
        bx = seg_arr[i, 1, 0]; by = seg_arr[i, 1, 1]
        if (abs(ax - px) < EPS and abs(ay - py) < EPS) or \
           (abs(bx - px) < EPS and abs(by - py) < EPS) or \
           (abs(ax - qx) < EPS and abs(ay - qy) < EPS) or \
           (abs(bx - qx) < EPS and abs(by - qy) < EPS):
            continue
        d1 = (qx - px) * (ay - py) - (qy - py) * (ax - px)
        d2 = (qx - px) * (by - py) - (qy - py) * (bx - px)
        d3 = (bx - ax) * (py - ay) - (by - ay) * (px - ax)
        d4 = (bx - ax) * (qy - ay) - (by - ay) * (qx - ax)
        if d1 * d2 < 0.0 and d3 * d4 < 0.0:
            return True
    return False


@njit(cache=True)
def _min_incident_angle(nx_, ny_, vx, vy, seg_arr):
    """Min angle (radians) between vector (vx,vy) and edges incident to node."""
    nv_norm = np.sqrt(vx * vx + vy * vy)
    if nv_norm < EPS:
        return np.pi
    best = np.pi
    for i in range(seg_arr.shape[0]):
        ax = seg_arr[i, 0, 0]; ay = seg_arr[i, 0, 1]
        bx = seg_arr[i, 1, 0]; by = seg_arr[i, 1, 1]
        if abs(ax - nx_) < EPS and abs(ay - ny_) < EPS:
            ox, oy = bx, by
        elif abs(bx - nx_) < EPS and abs(by - ny_) < EPS:
            ox, oy = ax, ay
        else:
            continue
        evx = ox - nx_; evy = oy - ny_
        ev_norm = np.sqrt(evx * evx + evy * evy)
        if ev_norm < EPS:
            continue
        cos_a = (vx * evx + vy * evy) / (nv_norm * ev_norm)
        if cos_a > 1.0:
            cos_a = 1.0
        elif cos_a < -1.0:
            cos_a = -1.0
        ang = np.arccos(cos_a)
        if ang < best:
            best = ang
    return best


def _connection_ok(mx, my, tx, ty, seg_arr, model, min_angle):
    """Visibility + (model 2) angle check for proposed edge (mx,my)→(tx,ty)."""
    if _segment_blocks(mx, my, tx, ty, seg_arr):
        return False
    if model == 2:
        vx, vy = tx - mx, ty - my
        if _min_incident_angle(mx, my, vx, vy, seg_arr) < min_angle:
            return False
        if _min_incident_angle(tx, ty, -vx, -vy, seg_arr) < min_angle:
            return False
    return True


# --------------------------------------------------------------------------- #
# Initialisation                                                              #
# --------------------------------------------------------------------------- #

def _initial_segments():
    """Unit square + 4 lines through the centre → 8 segments."""
    c = (CENTER_X, CENTER_Y)
    spokes = [
        ((0.0, 0.5), (1.0, 0.5)),
        ((0.5, 0.0), (0.5, 1.0)),
        ((0.0, 0.0), (1.0, 1.0)),
        ((1.0, 0.0), (0.0, 1.0)),
    ]
    out = []
    for p1, p2 in spokes:
        out.append((p1, c))
        out.append((c, p2))
    return out


def nkey(p, decimals=8):
    return (round(p[0], decimals), round(p[1], decimals))


# --------------------------------------------------------------------------- #
# Main simulation                                                             #
# --------------------------------------------------------------------------- #

def run_model(model=1, min_length=MIN_LENGTH, max_iter=200_000, seed=42,
              min_angle=np.pi / 4, area_coeff=0.05, area_scale=1.0,
              verbose=True):
    """Run the space-filling model.

    Parameters
    ----------
    model      : 1 (basic) or 2 (area + angle biases).
    min_length : segments shorter than this are never split.
    max_iter   : hard iteration cap.
    seed       : RNG seed.
    min_angle  : (model 2) minimum allowed intersection angle, radians.
    area_coeff : (model 2) coefficient in area_coeff·exp(-area_scale/r).
    area_scale : (model 2) scale factor in the area-constraint exponent.
    verbose    : print convergence / cap messages.

    Returns a list of segments [((x1,y1),(x2,y2)), ...].
    """
    rng = np.random.default_rng(seed)
    init = _initial_segments()

    # Segments stored in a capacity-doubling numpy buffer; live count is `n`.
    capacity = max(len(init) * 4, 64)
    seg_arr = np.empty((capacity, 2, 2), dtype=np.float64)
    for i, (a, b) in enumerate(init):
        seg_arr[i, 0] = a
        seg_arr[i, 1] = b
    n = len(init)

    # Node table: nkey → row index in node_arr (capacity-doubled).
    node_keys = {}
    node_arr = np.empty((capacity * 2, 2), dtype=np.float64)
    n_nodes = 0

    def add_node(pt):
        nonlocal n_nodes, node_arr
        k = nkey(pt)
        idx = node_keys.get(k)
        if idx is not None:
            return idx
        if n_nodes >= node_arr.shape[0]:
            new = np.empty((node_arr.shape[0] * 2, 2), dtype=np.float64)
            new[:n_nodes] = node_arr[:n_nodes]
            node_arr = new
        node_arr[n_nodes] = pt
        node_keys[k] = n_nodes
        n_nodes += 1
        return n_nodes - 1

    def grow_segs(extra):
        nonlocal seg_arr
        if n + extra <= seg_arr.shape[0]:
            return
        new_cap = max(seg_arr.shape[0] * 2, n + extra)
        new = np.empty((new_cap, 2, 2), dtype=np.float64)
        new[:n] = seg_arr[:n]
        seg_arr = new

    for a, b in init:
        add_node(a)
        add_node(b)

    converged = False
    for it in range(max_iter):
        sa = seg_arr[:n]
        dx = sa[:, 1, 0] - sa[:, 0, 0]
        dy = sa[:, 1, 1] - sa[:, 0, 1]
        lengths = np.hypot(dx, dy)

        eligible = lengths > min_length
        if model == 2:
            mx = (sa[:, 0, 0] + sa[:, 1, 0]) * 0.5
            my = (sa[:, 0, 1] + sa[:, 1, 1]) * 0.5
            r = np.hypot(mx - CENTER_X, my - CENTER_Y)
            with np.errstate(divide='ignore', invalid='ignore'):
                min_area = np.where(r > EPS,
                                    area_coeff * np.exp(-area_scale / r), 0.0)
            eligible &= (lengths * 0.5) ** 2 >= min_area

        elig_idx = np.flatnonzero(eligible)
        if elig_idx.size == 0:
            converged = True
            break

        # Length-weighted sampling via inverse-CDF (faster than np.random.choice).
        elig_lengths = lengths[elig_idx]
        cdf = np.cumsum(elig_lengths)
        seg_i = int(elig_idx[np.searchsorted(cdf, rng.random() * cdf[-1])])

        ax, ay = sa[seg_i, 0]
        bx, by = sa[seg_i, 1]
        midx, midy = 0.5 * (ax + bx), 0.5 * (ay + by)

        # Replace segment with its two halves (swap-remove + append).
        grow_segs(2)
        seg_arr[seg_i] = seg_arr[n - 1]
        n -= 1
        seg_arr[n, 0] = (ax, ay); seg_arr[n, 1] = (midx, midy); n += 1
        seg_arr[n, 0] = (midx, midy); seg_arr[n, 1] = (bx, by); n += 1

        m_idx = add_node((midx, midy))
        a_idx = node_keys[nkey((ax, ay))]
        b_idx = node_keys[nkey((bx, by))]

        # Candidate nodes for the new connection: all except m, a, b.
        keep = np.ones(n_nodes, dtype=bool)
        keep[[m_idx, a_idx, b_idx]] = False
        candidates = node_arr[:n_nodes][keep]
        if len(candidates) == 0:
            continue

        sa = seg_arr[:n]
        tree = cKDTree(candidates)
        k = min(30, len(candidates))
        _, idxs = tree.query((midx, midy), k=k)
        idxs = np.atleast_1d(idxs)

        best = _find_best(midx, midy, candidates, idxs, sa, model, min_angle)
        if best is None and k < len(candidates):
            _, idxs_all = tree.query((midx, midy), k=len(candidates))
            best = _find_best(midx, midy, candidates, idxs_all[k:], sa,
                              model, min_angle)

        if best is not None:
            grow_segs(1)
            seg_arr[n, 0] = (midx, midy)
            seg_arr[n, 1] = best
            n += 1
            add_node(best)

    if verbose:
        if converged:
            print(f"  Model {model}: converged at iteration {it} with {n} segments")
        else:
            print(f"  Model {model}: hit max_iter={max_iter}, {n} segments")

    return [(tuple(seg_arr[i, 0]), tuple(seg_arr[i, 1])) for i in range(n)]


def _find_best(mx, my, candidates, idxs, seg_arr, model, min_angle):
    """Return first candidate (by nearest-first order) that passes checks."""
    for idx in idxs:
        tx, ty = candidates[idx]
        if _connection_ok(mx, my, tx, ty, seg_arr, model, min_angle):
            return (float(tx), float(ty))
    return None


# --------------------------------------------------------------------------- #
# Parallel ensemble                                                           #
# --------------------------------------------------------------------------- #

def _ensemble_worker(args):
    kwargs, seed = args
    return run_model(seed=seed, verbose=False, **kwargs)


def run_ensemble(seeds, processes=None, **kwargs):
    """Run `run_model` over many seeds in parallel.

    Parameters
    ----------
    seeds     : iterable of RNG seeds.
    processes : worker count (default: os.cpu_count()).
    **kwargs  : forwarded to `run_model` (model, min_length, max_iter, ...).

    Returns a list of segment lists, one per seed (in input order).
    """
    seeds = list(seeds)
    if processes is None:
        processes = os.cpu_count() or 1
    payload = [(kwargs, s) for s in seeds]
    with Pool(processes=processes) as pool:
        return pool.map(_ensemble_worker, payload)


# --------------------------------------------------------------------------- #
# GeoDataFrame / shapefile export                                              #
# --------------------------------------------------------------------------- #

def compute_degrees(segments):
    deg = defaultdict(int)
    for a, b in segments:
        deg[nkey(a)] += 1
        deg[nkey(b)] += 1
    return dict(deg)


def to_geodataframe(segments, model_id):
    """Convert segments to a GeoDataFrame with length and node degrees."""
    degrees = compute_degrees(segments)
    records = [
        {
            'geometry':  LineString([a, b]),
            'model':     model_id,
            'length':    float(np.hypot(b[0] - a[0], b[1] - a[1])),
            'deg_start': degrees.get(nkey(a), 0),
            'deg_end':   degrees.get(nkey(b), 0),
        }
        for a, b in segments
    ]
    return gpd.GeoDataFrame(records)


def save_shapefile(segs, model_id, path):
    gdf = to_geodataframe(segs, model_id)
    gdf.to_file(path)
    print(f"  Saved {len(gdf)} segments → {path}")
    return gdf


def get_geodataframe(model=1, **kwargs):
    """Run `run_model` and return the result as a GeoDataFrame."""
    segs = run_model(model=model, **kwargs)
    return to_geodataframe(segs, model_id=model)


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #

if __name__ == '__main__':
    print("=== Space-Filling Street Network Model ===")
    print(f"  numba JIT: {'on' if HAS_NUMBA else 'off (install numba for speedup)'}\n")

    print("Running Model 1 (basic)...")
    segs1 = run_model(model=1, seed=42)
    save_shapefile(segs1, model_id=1, path='model1.shp')

    print("Running Model 2 (with area & angle biases)...")
    segs2 = run_model(model=2, seed=42)
    save_shapefile(segs2, model_id=2, path='model2.shp')
