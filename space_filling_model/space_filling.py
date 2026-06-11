#!/usr/bin/env python3
"""
Space-filling street network model.

Model 1 (basic): pick a segment with probability proportional to length, split
  at midpoint, connect midpoint to the closest visible intersection.

Model 2 (biased): same as Model 1 but only splits segments whose half-length²
  exceeds area_coeff·exp(-area_scale/r) (larger blocks near the periphery), and
  rejects connections whose nearest-edge angle falls outside
  [min_angle, max_angle] (max_angle=None disables the upper bound).

Both start from a unit square divided by 4 lines (horizontal, vertical, two
diagonals) crossing at the centre.

Segment selection is configurable via ``select`` ('length' default, 'uniform',
'power' with ``select_power``, or 'longest'); see ``run_model``.

Set ``p_stub`` (default 0.0) to inject random dead-end "stub" streets each
iteration via a separate RNG stream (so p_stub=0.0 is byte-identical to a run
without it). ``stub_len_max``/``snap_tol`` tune stub length and the T-in
distance; ``stub_at_midpoint_prob`` (default 0.5) splits stub anchoring between
the fresh midpoint and a random degree>=3 junction; ``stub_splittable`` (default
True) makes stubs first-class streets that feed back into growth, or inert
dead-ends when False.

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


def _nearest_segment(fx, fy, seg_arr, max_dist, exclude_xy=None):
    """Nearest segment to point (fx,fy) within `max_dist`.

    Returns (idx, foot_x, foot_y): the index into `seg_arr` of the closest
    segment and the foot of the perpendicular from the point onto it. If no
    segment lies within `max_dist`, returns (-1, fx, fy). Segments sharing the
    `exclude_xy` endpoint (the stub anchor) are never selected. (Local copy of
    the real_space_filling helper to avoid a circular import.)
    """
    if seg_arr.shape[0] == 0:
        return -1, fx, fy
    ax = seg_arr[:, 0, 0]; ay = seg_arr[:, 0, 1]
    bx = seg_arr[:, 1, 0]; by = seg_arr[:, 1, 1]
    dx = bx - ax; dy = by - ay
    L2 = dx * dx + dy * dy
    L2 = np.where(L2 < EPS, EPS, L2)
    t = ((fx - ax) * dx + (fy - ay) * dy) / L2
    t = np.clip(t, 0.0, 1.0)
    projx = ax + t * dx; projy = ay + t * dy
    dist = np.hypot(projx - fx, projy - fy)
    if exclude_xy is not None:
        ex, ey = exclude_xy
        share = (((np.abs(ax - ex) < EPS) & (np.abs(ay - ey) < EPS)) |
                 ((np.abs(bx - ex) < EPS) & (np.abs(by - ey) < EPS)))
        dist = np.where(share, np.inf, dist)
    i = int(np.argmin(dist))
    if dist[i] <= max_dist:
        return i, float(projx[i]), float(projy[i])
    return -1, fx, fy


def _connection_ok(mx, my, tx, ty, seg_arr, model, min_angle, max_angle=None):
    """Visibility + (model 2) angle check for proposed edge (mx,my)→(tx,ty).

    Model 2 rejects connections whose nearest-edge angle falls outside
    [min_angle, max_angle]; max_angle=None disables the upper bound.
    """
    if _segment_blocks(mx, my, tx, ty, seg_arr):
        return False
    if model == 2:
        vx, vy = tx - mx, ty - my
        a_m = _min_incident_angle(mx, my, vx, vy, seg_arr)
        if a_m < min_angle or (max_angle is not None and a_m > max_angle):
            return False
        a_t = _min_incident_angle(tx, ty, -vx, -vy, seg_arr)
        if a_t < min_angle or (max_angle is not None and a_t > max_angle):
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
              min_angle=np.pi / 4, max_angle=None, area_coeff=0.05,
              area_scale=1.0, select='length', select_power=1.0,
              p_stub=0.0, stub_len_max=None, snap_tol=None,
              stub_splittable=True, stub_at_midpoint_prob=0.5,
              noise_frac=0.0, noise_len_max=None, noise_snap_tol=None,
              noise_branch_prob=0.3, verbose=True):
    """Run the space-filling model.

    Parameters
    ----------
    model      : 1 (basic) or 2 (area + angle biases).
    min_length : segments shorter than this are never split.
    max_iter   : hard iteration cap.
    seed       : RNG seed.
    min_angle  : (model 2) minimum allowed intersection angle, radians.
    max_angle  : (model 2) maximum allowed nearest-edge angle, radians; None
                 disables the upper bound.
    area_coeff : (model 2) coefficient in area_coeff·exp(-area_scale/r).
    area_scale : (model 2) scale factor in the area-constraint exponent.
    select     : segment-selection strategy. 'length' (default) picks with
                 probability proportional to length; 'uniform' picks any
                 eligible segment equally; 'power' uses length**select_power;
                 'longest' deterministically splits the longest eligible segment.
    select_power : exponent for select='power' (1.0 ≡ 'length', 0.0 ≡ 'uniform').
    p_stub     : per-iteration probability of additionally emitting a random
                 dead-end "stub" street (0.0 = current behavior). Stubs use a
                 separate RNG stream, so p_stub=0.0 is byte-identical to a run
                 without it. Each stub grows from the fresh split midpoint or a
                 random degree>=3 intersection at a random heading/length and
                 either dangles as a degree-1 cul-de-sac or T's into a nearby
                 street within `snap_tol` (degree-2). Stub endpoints are confined
                 to the unit square and stubs never cross existing streets.
    stub_len_max : upper bound of stub length (default 3*min_length).
    snap_tol   : distance under which a stub free end T's into a nearby street
                 (default min_length). A small value yields mostly cul-de-sacs.
    stub_splittable : if True (default) stubs are normal streets — splittable and
                 valid connection targets — so they feed back into the
                 generative process. If False, stubs are non-splittable dead-ends
                 and dangling tips are never connection targets.
    stub_at_midpoint_prob : probability a stub is anchored at the fresh split
                 midpoint (default 0.5); otherwise it is anchored at a random
                 existing degree>=3 intersection. Falls back to the midpoint when
                 no degree>=3 node exists yet. 1.0 = always the new midpoint,
                 0.0 = always an existing junction (when one exists).
    verbose    : print convergence / cap messages.

    Returns a list of segments [((x1,y1),(x2,y2)), ...].

    Raises
    ------
    ValueError
        If `p_stub` or `stub_at_midpoint_prob` is outside [0, 1], or `select`
        is not a recognised strategy.
    """
    if not 0.0 <= p_stub <= 1.0:
        raise ValueError(f"p_stub must be in [0, 1], got {p_stub!r}")
    if not 0.0 <= stub_at_midpoint_prob <= 1.0:
        raise ValueError(
            f"stub_at_midpoint_prob must be in [0, 1], got "
            f"{stub_at_midpoint_prob!r}")
    if noise_frac < 0.0:
        raise ValueError(f"noise_frac must be >= 0, got {noise_frac!r}")
    if not 0.0 <= noise_branch_prob <= 1.0:
        raise ValueError(
            f"noise_branch_prob must be in [0, 1], got {noise_branch_prob!r}")

    rng = np.random.default_rng(seed)
    # Separate RNG for stubs so the main stream (segment selection and normal
    # connections) is identical regardless of p_stub; stubs are then a
    # reproducible additive layer. Placed stubs DO block later connections
    # (the network stays planar).
    stub_rng = np.random.default_rng(seed ^ 0xDEADBEEF) if p_stub > 0.0 else None
    if stub_len_max is None:
        stub_len_max = 3.0 * min_length
    if snap_tol is None:
        snap_tol = min_length
    init = _initial_segments()

    # Segments stored in a capacity-doubling numpy buffer; live count is `n`.
    # `splittable` marks segments eligible for midpoint splitting (stubs may be
    # added non-splittable).
    capacity = max(len(init) * 4, 64)
    seg_arr = np.empty((capacity, 2, 2), dtype=np.float64)
    splittable = np.zeros(capacity, dtype=bool)
    for i, (a, b) in enumerate(init):
        seg_arr[i, 0] = a
        seg_arr[i, 1] = b
        splittable[i] = True
    n = len(init)

    # Node table: nkey → row index in node_arr (capacity-doubled). `node_deg`
    # tracks primal degree (for degree>=3 stub anchors); `stub_node_set` holds
    # dangling stub tips that must not become connection targets.
    node_keys = {}
    node_arr = np.empty((capacity * 2, 2), dtype=np.float64)
    node_deg = np.zeros(capacity * 2, dtype=np.int64)
    n_nodes = 0
    deg3_idx = []
    stub_node_set = set()

    def add_node(pt):
        nonlocal n_nodes, node_arr, node_deg
        k = nkey(pt)
        idx = node_keys.get(k)
        if idx is not None:
            return idx
        if n_nodes >= node_arr.shape[0]:
            new = np.empty((node_arr.shape[0] * 2, 2), dtype=np.float64)
            new[:n_nodes] = node_arr[:n_nodes]
            node_arr = new
            new_deg = np.zeros(node_deg.shape[0] * 2, dtype=np.int64)
            new_deg[:n_nodes] = node_deg[:n_nodes]
            node_deg = new_deg
        node_arr[n_nodes] = pt
        node_keys[k] = n_nodes
        n_nodes += 1
        return n_nodes - 1

    def _bump(idx, amt):
        # Degrees only ever increase, so a node enters deg3_idx exactly once,
        # when it first reaches degree 3.
        old = int(node_deg[idx])
        node_deg[idx] = old + amt
        if old < 3 <= old + amt:
            deg3_idx.append(int(idx))

    def grow_segs(extra):
        nonlocal seg_arr, splittable
        if n + extra <= seg_arr.shape[0]:
            return
        new_cap = max(seg_arr.shape[0] * 2, n + extra)
        new = np.empty((new_cap, 2, 2), dtype=np.float64)
        new[:n] = seg_arr[:n]
        seg_arr = new
        new_s = np.zeros(new_cap, dtype=bool)
        new_s[:n] = splittable[:n]
        splittable = new_s

    for a, b in init:
        ia = add_node(a); ib = add_node(b)
        _bump(ia, 1); _bump(ib, 1)

    if select not in ('length', 'uniform', 'power', 'longest'):
        raise ValueError(f"unknown select={select!r}")

    converged = False
    for it in range(max_iter):
        sa = seg_arr[:n]
        dx = sa[:, 1, 0] - sa[:, 0, 0]
        dy = sa[:, 1, 1] - sa[:, 0, 1]
        lengths = np.hypot(dx, dy)

        eligible = splittable[:n] & (lengths > min_length)
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

        elig_lengths = lengths[elig_idx]
        seg_i = _select_segment(rng, elig_idx, elig_lengths, select,
                                select_power)

        ax, ay = sa[seg_i, 0]
        bx, by = sa[seg_i, 1]
        midx, midy = 0.5 * (ax + bx), 0.5 * (ay + by)

        # Replace segment with its two halves (swap-remove + append). Both
        # halves are splittable; the midpoint now carries them (degree += 2).
        grow_segs(2)
        last_split = splittable[n - 1]
        seg_arr[seg_i] = seg_arr[n - 1]
        splittable[seg_i] = last_split
        n -= 1
        seg_arr[n, 0] = (ax, ay); seg_arr[n, 1] = (midx, midy)
        splittable[n] = True; n += 1
        seg_arr[n, 0] = (midx, midy); seg_arr[n, 1] = (bx, by)
        splittable[n] = True; n += 1

        m_idx = add_node((midx, midy))
        a_idx = node_keys[nkey((ax, ay))]
        b_idx = node_keys[nkey((bx, by))]
        _bump(m_idx, 2)

        # --- additive stub: occasionally also emit a random dead-end stub from
        # this iteration's midpoint (or a random degree>=3 intersection). The
        # normal connection below still runs; placed stubs stay planar (they
        # block later connections). ---
        if p_stub > 0.0 and stub_rng.random() < p_stub:
            use_mid = (len(deg3_idx) == 0) or \
                (stub_rng.random() < stub_at_midpoint_prob)
            anchor_idx = m_idx if use_mid else \
                deg3_idx[int(stub_rng.integers(len(deg3_idx)))]
            ax0, ay0 = node_arr[anchor_idx]
            theta = stub_rng.random() * 2.0 * np.pi
            length = min_length + stub_rng.random() * (stub_len_max - min_length)
            fx = ax0 + length * np.cos(theta)
            fy = ay0 + length * np.sin(theta)

            # Optional T-in: snap the free end onto a nearby street (degree-2),
            # excluding edges incident to the anchor.
            snap_i, fx, fy = _nearest_segment(fx, fy, seg_arr[:n], snap_tol,
                                              exclude_xy=(ax0, ay0))

            ok = (0.0 <= fx <= 1.0) and (0.0 <= fy <= 1.0)
            if ok and _segment_blocks(ax0, ay0, fx, fy, seg_arr[:n]):
                ok = False
            if ok and model == 2:
                vx, vy = fx - ax0, fy - ay0
                a_m = _min_incident_angle(ax0, ay0, vx, vy, seg_arr[:n])
                if a_m < min_angle or (max_angle is not None and a_m > max_angle):
                    ok = False

            if ok:
                if snap_i >= 0:
                    # Split the snapped street at the contact point so a real
                    # junction is created (stub gains degree 2).
                    ux, uy = seg_arr[snap_i, 0]
                    wx, wy = seg_arr[snap_i, 1]
                    s_split = splittable[snap_i]
                    grow_segs(2)
                    seg_arr[snap_i] = seg_arr[n - 1]
                    splittable[snap_i] = splittable[n - 1]
                    n -= 1
                    c_idx = add_node((fx, fy))
                    seg_arr[n, 0] = (ux, uy); seg_arr[n, 1] = (fx, fy)
                    splittable[n] = s_split; n += 1
                    seg_arr[n, 0] = (fx, fy); seg_arr[n, 1] = (wx, wy)
                    splittable[n] = s_split; n += 1
                    _bump(c_idx, 2)
                    far_idx = c_idx
                else:
                    far_idx = add_node((fx, fy))

                grow_segs(1)
                seg_arr[n, 0] = (ax0, ay0); seg_arr[n, 1] = (fx, fy)
                splittable[n] = stub_splittable; n += 1
                _bump(anchor_idx, 1)
                _bump(far_idx, 1)
                if snap_i < 0 and not stub_splittable:
                    stub_node_set.add(far_idx)  # keep the cul-de-sac dangling

        # Candidate nodes for the new connection: all except m, a, b (and any
        # excluded dangling stub tips).
        keep = np.ones(n_nodes, dtype=bool)
        keep[[m_idx, a_idx, b_idx]] = False
        if stub_node_set:
            keep[list(stub_node_set)] = False
        candidates = node_arr[:n_nodes][keep]
        if len(candidates) == 0:
            continue

        sa = seg_arr[:n]
        tree = cKDTree(candidates)
        k = min(30, len(candidates))
        _, idxs = tree.query((midx, midy), k=k)
        idxs = np.atleast_1d(idxs)

        best = _find_best(midx, midy, candidates, idxs, sa, model, min_angle,
                          max_angle)
        if best is None and k < len(candidates):
            _, idxs_all = tree.query((midx, midy), k=len(candidates))
            best = _find_best(midx, midy, candidates, idxs_all[k:], sa,
                              model, min_angle, max_angle)

        if best is not None:
            grow_segs(1)
            seg_arr[n, 0] = (midx, midy)
            seg_arr[n, 1] = best
            splittable[n] = True
            n += 1
            best_idx = add_node(best)
            _bump(m_idx, 1)        # midpoint now a T-junction (degree 3)
            _bump(best_idx, 1)

    if verbose:
        if converged:
            print(f"  Model {model}: converged at iteration {it} with {n} segments")
        else:
            print(f"  Model {model}: hit max_iter={max_iter}, {n} segments")

    segments = [(tuple(seg_arr[i, 0]), tuple(seg_arr[i, 1])) for i in range(n)]
    if noise_frac > 0.0:
        # Post-growth decoration on its own RNG stream: noise_frac=0.0 is
        # byte-identical to a run without it.
        segments = add_noise_streets(
            segments, noise_frac=noise_frac, seed=seed ^ 0xBADC0FFE,
            min_length=min_length, noise_len_max=noise_len_max,
            snap_tol=noise_snap_tol, branch_prob=noise_branch_prob,
            min_angle=min_angle, verbose=verbose)
    return segments


def _find_best(mx, my, candidates, idxs, seg_arr, model, min_angle,
               max_angle=None):
    """Return first candidate (by nearest-first order) that passes checks."""
    for idx in idxs:
        tx, ty = candidates[idx]
        if _connection_ok(mx, my, tx, ty, seg_arr, model, min_angle, max_angle):
            return (float(tx), float(ty))
    return None


def _select_segment(rng, elig_idx, elig_lengths, select='length',
                    select_power=1.0):
    """Return the chosen segment index (an int drawn from elig_idx).

    Strategies
    ----------
    'length'  : probability proportional to length (default; current behaviour).
    'uniform' : every eligible segment equally likely.
    'power'   : probability proportional to length ** select_power.
    'longest' : deterministic argmax (ties -> lowest elig_idx).
    """
    if select == 'length':
        cdf = np.cumsum(elig_lengths)
        return int(elig_idx[np.searchsorted(cdf, rng.random() * cdf[-1])])
    if select == 'uniform':
        return int(rng.choice(elig_idx))
    if select == 'power':
        cdf = np.cumsum(elig_lengths ** select_power)
        return int(elig_idx[np.searchsorted(cdf, rng.random() * cdf[-1])])
    if select == 'longest':
        return int(elig_idx[np.argmax(elig_lengths)])
    raise ValueError(f"unknown select={select!r}")


# --------------------------------------------------------------------------- #
# Post-growth low-degree noise decoration                                     #
# --------------------------------------------------------------------------- #

def add_noise_streets(segments, noise_frac=0.5, seed=0,
                      min_length=MIN_LENGTH, noise_len_max=None,
                      snap_tol=None, branch_prob=0.3,
                      min_angle=np.pi / 4, max_tries=50, verbose=False):
    """Decorate a finished network with random low-degree "noise" streets.

    Post-growth counterpart to ``p_stub``: takes the segment list returned by
    `run_model` and adds ``round(noise_frac * len(segments))`` short streets,
    each anchored mid-block on an existing street (the host is split at the
    anchor, so DGDC sees the noise street touch exactly one stroke there).
    Because growth never builds on them, noise streets keep low dual-graph
    degrees: a dangling free end gives degree 1, a snapped (T-in) end gives
    degree 2, and each child anchored on a noise street (``branch_prob``)
    adds +1 to its parent — so the degree-1/2/3 population is directly
    tunable. Sweep ``noise_frac`` on one saved backbone without re-running
    growth.

    Parameters
    ----------
    segments    : list of ((x1,y1),(x2,y2)) street segments (unit square).
    noise_frac  : number of noise streets to add as a fraction of
                  ``len(segments)``; unbounded above, 0.0 is a no-op.
    seed        : RNG seed of the decoration stream.
    min_length  : lower bound of noise street length.
    noise_len_max : upper bound of noise street length (default 3*min_length).
    snap_tol    : free ends within this distance of a street T into it
                  (splitting it, dual degree 2); 0 disables snapping so all
                  noise streets dangle (degree 1). Default min_length.
    branch_prob : probability of anchoring on a previously added noise street
                  (small trees -> parents reach degree 2-3). Falls back to the
                  general pool while no noise street exists yet.
    min_angle   : minimum angle (radians) between a noise street and its
                  host's axis; keep above the DGDC a_threshold so noise is
                  never merged into the host stroke.
    max_tries   : placement attempts per noise street before it is skipped.
    verbose     : print a placement summary.

    Returns a new segment list: the (re-split) input plus the noise streets.
    The result is planar (noise never crosses existing streets) and confined
    to the unit square.

    Raises
    ------
    ValueError
        If ``noise_frac`` is negative or ``branch_prob`` is outside [0, 1].
    """
    if noise_frac < 0.0:
        raise ValueError(f"noise_frac must be >= 0, got {noise_frac!r}")
    if not 0.0 <= branch_prob <= 1.0:
        raise ValueError(f"branch_prob must be in [0, 1], got {branch_prob!r}")
    if noise_len_max is None:
        noise_len_max = 3.0 * min_length
    if snap_tol is None:
        snap_tol = min_length

    target = int(round(noise_frac * len(segments)))
    if target == 0:
        return list(segments)

    rng = np.random.default_rng(seed)

    n = len(segments)
    capacity = max(n + 3 * target, 64)
    seg_arr = np.empty((capacity, 2, 2), dtype=np.float64)
    is_noise = np.zeros(capacity, dtype=bool)
    for i, (a, b) in enumerate(segments):
        seg_arr[i, 0] = a
        seg_arr[i, 1] = b

    def grow(extra):
        nonlocal seg_arr, is_noise
        if n + extra <= seg_arr.shape[0]:
            return
        new_cap = max(seg_arr.shape[0] * 2, n + extra)
        new = np.empty((new_cap, 2, 2), dtype=np.float64)
        new[:n] = seg_arr[:n]
        seg_arr = new
        new_m = np.zeros(new_cap, dtype=bool)
        new_m[:n] = is_noise[:n]
        is_noise = new_m

    def split_row(j, px, py):
        # Overwrite row j with (a, p) and append (p, b): all other row indices
        # stay stable. Both halves inherit the row's noise flag.
        nonlocal n
        bx_, by_ = seg_arr[j, 1]
        grow(1)
        seg_arr[j, 1] = (px, py)
        seg_arr[n, 0] = (px, py)
        seg_arr[n, 1] = (bx_, by_)
        is_noise[n] = is_noise[j]
        n += 1

    # Hiding a row behind this far-away degenerate sentinel excludes it from
    # _segment_blocks / _nearest_segment in O(1): a point lying *on* a row
    # (anchor on host, snapped foot on target) makes the proper-crossing test
    # numerically unstable, so those rows must not take part in the checks.
    _far = ((-10.0, -10.0), (-10.0, -10.0))
    GUARD = 1e-6  # degeneracy guard, safely above nkey's 1e-8 rounding

    placed = 0
    skipped = 0
    for _ in range(target):
        success = False
        for _try in range(max_tries):
            # ---- parent street (length-weighted; optionally noise-only) ----
            use_noise_pool = rng.random() < branch_prob
            if use_noise_pool and is_noise[:n].any():
                pool = np.flatnonzero(is_noise[:n])
            else:
                pool = np.arange(n)
            sa = seg_arr[:n]
            lens = np.hypot(sa[pool, 1, 0] - sa[pool, 0, 0],
                            sa[pool, 1, 1] - sa[pool, 0, 1])
            cdf = np.cumsum(lens)
            h = int(pool[np.searchsorted(cdf, rng.random() * cdf[-1])])
            hax, hay = seg_arr[h, 0]
            hbx, hby = seg_arr[h, 1]
            if np.hypot(hbx - hax, hby - hay) < 10.0 * GUARD:
                continue  # host too short to split cleanly
            host_phi = np.arctan2(hby - hay, hbx - hax)

            # ---- anchor, heading, length ----
            t = 0.1 + 0.8 * rng.random()
            ax0 = hax + t * (hbx - hax)
            ay0 = hay + t * (hby - hay)
            theta = rng.random() * 2.0 * np.pi
            d = (theta - host_phi) % np.pi
            if min(d, np.pi - d) < min_angle:
                continue  # too close to the host axis (DGDC would merge)
            length = min_length + rng.random() * (noise_len_max - min_length)
            fx = ax0 + length * np.cos(theta)
            fy = ay0 + length * np.sin(theta)

            # ---- snap & validity (host, and snapped row, hidden) ----
            saved_h = seg_arr[h].copy()
            seg_arr[h] = _far
            snap_j = -1
            endpoint_connect = False
            saved_j = None
            if snap_tol > 0.0:
                snap_j, fx, fy = _nearest_segment(fx, fy, seg_arr[:n],
                                                  snap_tol)
            if snap_j >= 0:
                ux, uy = seg_arr[snap_j, 0]
                wx, wy = seg_arr[snap_j, 1]
                if np.hypot(fx - ux, fy - uy) < GUARD:
                    fx, fy = float(ux), float(uy)   # join the existing node
                    endpoint_connect = True
                elif np.hypot(fx - wx, fy - wy) < GUARD:
                    fx, fy = float(wx), float(wy)
                    endpoint_connect = True
                else:
                    saved_j = seg_arr[snap_j].copy()
                    seg_arr[snap_j] = _far
            # Snapping moved the free end, so re-check the host-axis angle on
            # the final geometry (otherwise a snapped end can land collinear
            # with — even on — the host).
            df = (np.arctan2(fy - ay0, fx - ax0) - host_phi) % np.pi
            ok = (0.0 <= fx <= 1.0 and 0.0 <= fy <= 1.0
                  and np.hypot(fx - ax0, fy - ay0) >= GUARD
                  and min(df, np.pi - df) >= min_angle
                  and not _segment_blocks(ax0, ay0, fx, fy, seg_arr[:n]))
            # Also verify both split halves of the host stay planar.  The
            # original full host segment was fine, but a shorter half can cross
            # segments that lay "beyond" the anchor along the host's axis.
            # The anchor endpoint is shared with both halves, so _segment_blocks
            # treats segments incident to (ax0,ay0) as non-blocking — that is
            # correct (shared node, not a proper crossing).
            #
            # An additional degenerate case: when the host is collinear with a
            # co-incident sibling (same endpoint, same direction), a short split
            # half can *overlap* the sibling.  _segment_blocks misses this (all
            # cross products ≈ 0 for collinear points, so d1*d2 ≈ 0 and the < 0
            # test fails).  Guard against overlapping halves by checking that
            # neither half's midpoint lies on any existing segment.
            if ok:
                ok = (not _segment_blocks(hax, hay, ax0, ay0, seg_arr[:n])
                      and not _segment_blocks(ax0, ay0, hbx, hby, seg_arr[:n]))
            if ok:
                # Overlap guard: midpoint of each half must not lie on any
                # existing segment (collinear-overlap rejection).
                m1x = 0.5 * (hax + ax0); m1y = 0.5 * (hay + ay0)
                m2x = 0.5 * (ax0 + hbx); m2y = 0.5 * (ay0 + hby)
                _ov_tol = GUARD
                j1, _, _ = _nearest_segment(m1x, m1y, seg_arr[:n], _ov_tol)
                j2, _, _ = _nearest_segment(m2x, m2y, seg_arr[:n], _ov_tol)
                if j1 >= 0 or j2 >= 0:
                    ok = False
            seg_arr[h] = saved_h
            if saved_j is not None:
                seg_arr[snap_j] = saved_j
            if not ok:
                continue

            # ---- commit: split host (and snapped street), append noise ----
            split_row(h, ax0, ay0)
            if snap_j >= 0 and not endpoint_connect:
                split_row(snap_j, fx, fy)
            grow(1)
            seg_arr[n, 0] = (ax0, ay0)
            seg_arr[n, 1] = (fx, fy)
            is_noise[n] = True
            n += 1
            success = True
            break
        if success:
            placed += 1
        else:
            skipped += 1

    if verbose:
        print(f"  noise: placed {placed}/{target} noise streets "
              f"({skipped} skipped)")

    return [(tuple(seg_arr[i, 0]), tuple(seg_arr[i, 1])) for i in range(n)]


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