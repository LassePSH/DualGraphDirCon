#!/usr/bin/env python3
"""
Real-area space-filling street network model (circular domain).

Wraps `space_filling.run_model` so that the unit-square domain is replaced by a
circular area in local UTM metres. The circle's area matches the convex-hull
area of a cached OSM city graph (if a `city` is given) or an explicit
`area_m2`. This lets model parameters be compared in metres against real-city
sizes without any coastline/water confounds.

Init is "spokes": the circle boundary plus N lines from the centre to the
boundary (analogue of the unit-square + diagonals init).

Geometry units throughout the public API are metres. The returned GeoDataFrame
is always georeferenced: when a `city` is given its local UTM CRS is used,
otherwise the CRS defaults to EPSG:32632 (UTM zone 32N). Defaults:
`min_length=50` m, `min_angle=π/4`.
"""

from __future__ import annotations

import os
import sys
from typing import Optional

import geopandas as gpd
import numpy as np
import osmnx as ox
from scipy.spatial import cKDTree
from shapely.geometry import LineString, Point, Polygon

sys.path.insert(0, os.path.dirname(__file__))
from space_filling import (  # noqa: E402  reuse hot inner loops
    _find_best,
    _segment_blocks,
    _min_incident_angle,
    EPS,
)

# Default CRS for an abstract circle that has no real-world location.
DEFAULT_CRS = 'EPSG:32632'  # UTM zone 32N (metres)


def nkey(p, decimals=6):
    """Coordinate key — 6 decimals (sub-mm) is robust to float64 roundtrip."""
    return (round(float(p[0]), decimals), round(float(p[1]), decimals))

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
GRAPH_DIR = os.path.join(REPO_ROOT, 'data', 'city_graphs')


# --------------------------------------------------------------------------- #
# Area construction                                                           #
# --------------------------------------------------------------------------- #

def _city_safe_name(city: str) -> str:
    return city.replace('/', '_').replace(' ', '_')


def load_city_graph(city: str):
    """Load cached graphml and project to local UTM. Returns projected graph."""
    path = os.path.join(GRAPH_DIR, _city_safe_name(city) + '.graphml')
    if not os.path.exists(path):
        # also accept already-safe name
        path = os.path.join(GRAPH_DIR, city + '.graphml')
    if not os.path.exists(path):
        raise FileNotFoundError(f'No cached graph for {city!r} in {GRAPH_DIR}')
    G = ox.load_graphml(path)
    return ox.project_graph(G)


def circular_area_polygon(area_m2: float, center=(0.0, 0.0), resolution: int = 128) -> Polygon:
    """Circle with prescribed area (m^2) centred at `center`."""
    radius = float(np.sqrt(area_m2 / np.pi))
    return Point(center).buffer(radius, resolution=resolution)


def area_polygon_from_graph(G_proj) -> Polygon:
    """Convex hull of all node coordinates in the projected graph."""
    xs = np.fromiter((d['x'] for _, d in G_proj.nodes(data=True)), dtype=float)
    ys = np.fromiter((d['y'] for _, d in G_proj.nodes(data=True)), dtype=float)
    pts = gpd.GeoSeries(gpd.points_from_xy(xs, ys))
    hull = pts.union_all().convex_hull
    if not isinstance(hull, Polygon):
        raise RuntimeError('Convex hull is not a polygon (too few nodes?)')
    return hull


# --------------------------------------------------------------------------- #
# Blocking segments (area boundary)                                           #
# --------------------------------------------------------------------------- #

def _polygon_boundary_segments(poly):
    """Yield (a, b) tuples for every consecutive vertex pair in a polygon."""
    if poly.is_empty:
        return
    rings = [poly.exterior, *poly.interiors] if hasattr(poly, 'exterior') else []
    for ring in rings:
        coords = list(ring.coords)
        for i in range(len(coords) - 1):
            yield (coords[i][:2], coords[i + 1][:2])


def boundary_segments(area_poly: Polygon) -> list:
    """Non-splittable blocking segments for the outer area boundary."""
    return list(_polygon_boundary_segments(area_poly))


# --------------------------------------------------------------------------- #
# Initialisation                                                              #
# --------------------------------------------------------------------------- #

def init_spokes(area_poly: Polygon, n_spokes: int = 8) -> list:
    """Centre → boundary spokes."""
    c = area_poly.representative_point()
    cx, cy = c.x, c.y
    boundary = area_poly.exterior
    out = []
    for i in range(n_spokes):
        t = i / n_spokes
        p = boundary.interpolate(t, normalized=True)
        out.append(((cx, cy), (p.x, p.y)))
    return out


# --------------------------------------------------------------------------- #
# Main simulation (real-area variant)                                         #
# --------------------------------------------------------------------------- #

def _nearest_segment(fx, fy, seg_arr, max_dist, exclude_xy=None):
    """Nearest segment to point (fx,fy) within `max_dist`.

    Returns (idx, foot_x, foot_y): the index into `seg_arr` of the closest
    segment and the foot of the perpendicular from the point onto it. If no
    segment lies within `max_dist`, returns (-1, fx, fy). Segments sharing the
    `exclude_xy` endpoint (the stub anchor) are never selected.
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


def _run(seg_init, blocking_segs, area_poly, *,
         model=1, min_length=50.0, max_iter=200_000, seed=42,
         min_angle=np.pi / 4, area_coeff=2500.0, area_scale=None,
         p_stub=0.0, stub_len_max=None, snap_tol=None,
         kdtree_rebuild_every=256, progress_every=0, verbose=True):
    """Core loop. `seg_init` are splittable seeds; `blocking_segs` are not."""
    rng = np.random.default_rng(seed)
    # Separate RNG for stubs so the main RNG stream (segment selection and
    # normal connections) is identical regardless of p_stub; stubs are then a
    # reproducible additive layer. Placed stubs DO block later connections
    # (the network stays planar), so the resulting backbone can differ slightly
    # from the p_stub=0 case even though the RNG draws match.
    stub_rng = np.random.default_rng(seed ^ 0xDEADBEEF) if p_stub > 0.0 else None
    if stub_len_max is None:
        stub_len_max = 3.0 * min_length
    if snap_tol is None:
        snap_tol = min_length

    n_init = len(seg_init)
    n_block = len(blocking_segs)
    capacity = max((n_init + n_block) * 4, 256)
    seg_arr = np.empty((capacity, 2, 2), dtype=np.float64)
    splittable = np.zeros(capacity, dtype=bool)
    is_street = np.zeros(capacity, dtype=bool)

    for i, (a, b) in enumerate(seg_init):
        seg_arr[i, 0] = a; seg_arr[i, 1] = b
        splittable[i] = True; is_street[i] = True
    for j, (a, b) in enumerate(blocking_segs):
        k = n_init + j
        seg_arr[k, 0] = a; seg_arr[k, 1] = b
        splittable[k] = False; is_street[k] = False
    n = n_init + n_block

    node_keys = {}
    node_arr = np.empty((capacity * 2, 2), dtype=np.float64)
    node_deg = np.zeros(capacity * 2, dtype=np.int64)
    n_nodes = 0
    deg3_idx = []  # node indices whose primal degree has reached >= 3
    stub_node_set = set()  # far-end indices of placed stubs (excluded from normal candidates)

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
        # Degrees only ever increase in this model, so a node enters deg3_idx
        # exactly once, when it first reaches degree 3.
        old = int(node_deg[idx])
        node_deg[idx] = old + amt
        if old < 3 <= old + amt:
            deg3_idx.append(int(idx))

    def grow(extra):
        nonlocal seg_arr, splittable, is_street
        if n + extra <= seg_arr.shape[0]:
            return
        new_cap = max(seg_arr.shape[0] * 2, n + extra)
        new = np.empty((new_cap, 2, 2), dtype=np.float64)
        new[:n] = seg_arr[:n]
        seg_arr = new
        new_s = np.zeros(new_cap, dtype=bool)
        new_s[:n] = splittable[:n]
        splittable = new_s
        new_st = np.zeros(new_cap, dtype=bool)
        new_st[:n] = is_street[:n]
        is_street = new_st

    for a, b in seg_init:
        ia = add_node(a); ib = add_node(b); _bump(ia, 1); _bump(ib, 1)
    for a, b in blocking_segs:
        ia = add_node(a); ib = add_node(b); _bump(ia, 1); _bump(ib, 1)

    # area-bias centring (model 2): use mean of initial seed nodes as a proxy
    # for the network centre.
    if n_init > 0:
        seed_pts = np.array([list(p) for ab in seg_init for p in ab])
        cx, cy = float(seed_pts[:, 0].mean()), float(seed_pts[:, 1].mean())
    else:
        c = area_poly.centroid
        cx, cy = c.x, c.y
    if area_scale is None:
        # characteristic radius of the area in metres
        area_scale = float(np.sqrt(area_poly.area / np.pi))

    # Persistent KD-tree state. Rebuilding the tree from scratch every iteration
    # was O(N log N) per step → O(N^2 log N) overall and dominated runtime. We
    # instead rebuild only every `kdtree_rebuild_every` node insertions and
    # brute-force the handful of nodes added since the last rebuild.
    tree = None
    tree_n = 0

    converged = False
    last_it = 0
    show_progress = bool(progress_every) and verbose
    for it in range(max_iter):
        last_it = it
        # Lightweight progress: a single rewritten line every `progress_every`
        # iterations (the modulo test is the only per-iteration cost).
        if show_progress and it % progress_every == 0:
            print(f'\r  Real model {model}: iter {it}/{max_iter}, '
                  f'{n} segs', end='', flush=True)
        sa = seg_arr[:n]
        dx = sa[:, 1, 0] - sa[:, 0, 0]
        dy = sa[:, 1, 1] - sa[:, 0, 1]
        lengths = np.hypot(dx, dy)

        eligible = splittable[:n] & (lengths > min_length)
        if model == 2:
            mx = (sa[:, 0, 0] + sa[:, 1, 0]) * 0.5
            my = (sa[:, 0, 1] + sa[:, 1, 1]) * 0.5
            r = np.hypot(mx - cx, my - cy)
            with np.errstate(divide='ignore', invalid='ignore'):
                min_area = np.where(r > EPS,
                                    area_coeff * np.exp(-area_scale / r), 0.0)
            eligible &= (lengths * 0.5) ** 2 >= min_area

        elig_idx = np.flatnonzero(eligible)
        if elig_idx.size == 0:
            converged = True
            break

        elig_lengths = lengths[elig_idx]
        cdf = np.cumsum(elig_lengths)
        seg_i = int(elig_idx[np.searchsorted(cdf, rng.random() * cdf[-1])])

        ax, ay = sa[seg_i, 0]
        bx, by = sa[seg_i, 1]
        midx, midy = 0.5 * (ax + bx), 0.5 * (ay + by)

        # split: swap-remove the parent, append two halves (both splittable)
        grow(2)
        last_split = splittable[n - 1]
        last_street = is_street[n - 1]
        seg_arr[seg_i] = seg_arr[n - 1]
        splittable[seg_i] = last_split
        is_street[seg_i] = last_street
        n -= 1
        seg_arr[n, 0] = (ax, ay); seg_arr[n, 1] = (midx, midy)
        splittable[n] = True; is_street[n] = True; n += 1
        seg_arr[n, 0] = (midx, midy); seg_arr[n, 1] = (bx, by)
        splittable[n] = True; is_street[n] = True; n += 1

        m_idx = add_node((midx, midy))
        a_idx = node_keys[nkey((ax, ay))]
        b_idx = node_keys[nkey((bx, by))]
        _bump(m_idx, 2)  # midpoint now carries the two halves

        # --- branch-or-stub: occasionally emit a dead-end stub instead of the
        # usual nearest-visible connection. Stubs are streets but never split. ---
        if p_stub > 0.0 and stub_rng.random() < p_stub:
            use_mid = (len(deg3_idx) == 0) or (stub_rng.random() < 0.5)
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

            ok = area_poly.contains(Point(fx, fy))
            if ok and _segment_blocks(ax0, ay0, fx, fy, seg_arr[:n]):
                ok = False
            if ok and model == 2:
                vx, vy = fx - ax0, fy - ay0
                if _min_incident_angle(ax0, ay0, vx, vy, seg_arr[:n]) < min_angle:
                    ok = False

            if ok:
                if snap_i >= 0:
                    # Split the snapped street at the contact point so DGDC sees
                    # a real junction (stub stroke gains degree 2).
                    ux, uy = seg_arr[snap_i, 0]
                    wx, wy = seg_arr[snap_i, 1]
                    s_split = splittable[snap_i]
                    s_street = is_street[snap_i]
                    grow(2)
                    seg_arr[snap_i] = seg_arr[n - 1]
                    splittable[snap_i] = splittable[n - 1]
                    is_street[snap_i] = is_street[n - 1]
                    n -= 1
                    c_idx = add_node((fx, fy))
                    seg_arr[n, 0] = (ux, uy); seg_arr[n, 1] = (fx, fy)
                    splittable[n] = s_split; is_street[n] = s_street; n += 1
                    seg_arr[n, 0] = (fx, fy); seg_arr[n, 1] = (wx, wy)
                    splittable[n] = s_split; is_street[n] = s_street; n += 1
                    _bump(c_idx, 2)
                    far_idx = c_idx
                else:
                    far_idx = add_node((fx, fy))

                grow(1)
                seg_arr[n, 0] = (ax0, ay0); seg_arr[n, 1] = (fx, fy)
                splittable[n] = False; is_street[n] = True; n += 1
                _bump(anchor_idx, 1)
                _bump(far_idx, 1)
                if snap_i < 0:  # dangling tip: never let a normal connection
                    stub_node_set.add(far_idx)  # target it (keeps the cul-de-sac)

        if n_nodes <= 3:
            continue

        sa = seg_arr[:n]

        # (Re)build the tree only when enough new nodes have accumulated.
        if tree is None or (n_nodes - tree_n) >= kdtree_rebuild_every:
            tree = cKDTree(node_arr[:n_nodes].copy())
            tree_n = n_nodes

        exclude = (m_idx, a_idx, b_idx)

        def _ordered(limit):
            """Candidate coords nearest-first, excluding m/a/b nodes."""
            k = min(limit, tree_n)
            if k > 0:
                d, gi = tree.query((midx, midy), k=k)
                gi = np.atleast_1d(gi).astype(np.intp)
                d = np.atleast_1d(d)
            else:
                gi = np.empty(0, dtype=np.intp)
                d = np.empty(0)
            if n_nodes > tree_n:  # brute-force the nodes added since last build
                pend = np.arange(tree_n, n_nodes, dtype=np.intp)
                pv = node_arr[pend]
                pd = np.hypot(pv[:, 0] - midx, pv[:, 1] - midy)
                gi = np.concatenate([gi, pend])
                d = np.concatenate([d, pd])
            gi = gi[np.argsort(d, kind='stable')]
            mask = (gi != exclude[0]) & (gi != exclude[1]) & (gi != exclude[2])
            if stub_node_set:
                mask &= np.array([int(i) not in stub_node_set for i in gi], dtype=bool)
            return node_arr[gi[mask]]

        cand = _ordered(30)
        if len(cand) == 0:
            continue
        best = _find_best(midx, midy, cand, range(len(cand)), sa,
                          model, min_angle)
        if best is None and tree_n > 30:  # widen search to all nodes
            cand_all = _ordered(tree_n)
            best = _find_best(midx, midy, cand_all, range(len(cand_all)),
                              sa, model, min_angle)

        if best is not None:
            grow(1)
            seg_arr[n, 0] = (midx, midy)
            seg_arr[n, 1] = best
            splittable[n] = True; is_street[n] = True
            n += 1
            best_idx = add_node(best)
            _bump(m_idx, 1)        # midpoint now a T-junction (degree 3)
            _bump(best_idx, 1)

    if show_progress:
        print()  # finish the in-place progress line
    if verbose:
        tag = 'converged' if converged else f'hit max_iter={max_iter}'
        print(f'  Real model {model}: {tag} at iter {last_it} with {n} segs '
              f'({n_block} blocking, {n - n_block} streets)')

    # return only the street segments; boundary blockers stay separate
    return [(tuple(seg_arr[i, 0]), tuple(seg_arr[i, 1]))
            for i in range(n) if is_street[i]]


def run_circular_model(city: Optional[str] = None, *, area_m2: Optional[float] = None,
                       model: int = 1, min_length: float = 50.0,
                       min_angle: float = np.pi / 4, max_iter: int = 200_000,
                       seed: int = 42, n_spokes: int = 8,
                       area_coeff: float = 2500.0,
                       area_scale: Optional[float] = None,
                       p_stub: float = 0.0,
                       stub_len_max: Optional[float] = None,
                       snap_tol: Optional[float] = None,
                       kdtree_rebuild_every: int = 256,
                       progress_every: int = 1000,
                       crs=None, verbose: bool = True):
    """Run the space-filling model on a plain circular domain.

    The circle's area matches the convex-hull area of `city` (if given) or the
    explicit `area_m2`. Useful for comparing model parameters in metres against
    real-city parameters without confounds from coastline geometry.

    Returns
    -------
    dict with keys:
        segments      : list of ((x1,y1),(x2,y2)) in metres.
        area_polygon  : shapely Polygon (metres).
        crs           : CRS of the geometry (never None).
        init, model, area_m2.
    """
    if area_m2 is None:
        if city is None:
            raise ValueError('Provide either `city` or `area_m2`')
        G_proj = load_city_graph(city)
        area_m2 = float(area_polygon_from_graph(G_proj).area)
        if crs is None:
            crs = G_proj.graph['crs']

    if crs is None:
        crs = DEFAULT_CRS

    area = circular_area_polygon(area_m2)
    seeds = init_spokes(area, n_spokes=n_spokes)
    blockers = boundary_segments(area)

    segments = _run(seeds, blockers, area,
                    model=model, min_length=min_length, max_iter=max_iter,
                    seed=seed, min_angle=min_angle,
                    area_coeff=area_coeff, area_scale=area_scale,
                    p_stub=p_stub, stub_len_max=stub_len_max, snap_tol=snap_tol,
                    kdtree_rebuild_every=kdtree_rebuild_every,
                    progress_every=progress_every, verbose=verbose)

    return {
        'segments': segments,
        'area_polygon': area,
        'crs': crs,
        'init': 'spokes',
        'model': model,
        'area_m2': area_m2,
    }


# --------------------------------------------------------------------------- #
# GeoDataFrame export                                                         #
# --------------------------------------------------------------------------- #

def to_geodataframe(result: dict, model_id: Optional[int] = None) -> gpd.GeoDataFrame:
    """Build a GeoDataFrame of street segments in the projected (metres) CRS."""
    from collections import defaultdict
    segments = result['segments']
    deg = defaultdict(int)
    for a, b in segments:
        deg[nkey(a)] += 1
        deg[nkey(b)] += 1
    mid = model_id if model_id is not None else result.get('model', 0)
    geoms = [LineString([a, b]) for a, b in segments]
    return gpd.GeoDataFrame(
        {
            'model':     mid,
            'length':    [float(np.hypot(b[0] - a[0], b[1] - a[1]))
                          for a, b in segments],
            'deg_start': [deg[nkey(a)] for a, b in segments],
            'deg_end':   [deg[nkey(b)] for a, b in segments],
        },
        geometry=geoms,
        crs=result['crs'],
    )


# --------------------------------------------------------------------------- #
# Entry point                                                                 #
# --------------------------------------------------------------------------- #

if __name__ == '__main__':
    area_m2 = np.pi * 500.0 ** 2
    print(f'=== Circular space-filling model: area={area_m2:.0f} m^2 ===')
    for model in (1, 2):
        print(f'\n[model {model}]')
        res = run_circular_model(area_m2=area_m2, model=model,
                                 min_length=80.0, max_iter=20_000, seed=42)
        gdf = to_geodataframe(res, model_id=model)
        out = f'circular_m{model}.gpkg'
        gdf.to_file(out, driver='GPKG')
        print(f'  wrote {out} ({len(gdf)} segments, crs={gdf.crs})')
