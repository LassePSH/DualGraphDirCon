#!/usr/bin/env python3
"""
Real-area space-filling street network model.

Wraps `space_filling.run_model` so that the unit-square domain is replaced by
the convex hull of a cached OSM city graph (in local UTM metres), and any
proposed edge that would cross water, coastline, or the area boundary is
rejected.

Two init modes:
  * "spokes" — area boundary + N lines from the centroid to the boundary
    (analogue of the unit-square + diagonals init).
  * "roads"  — boundary + the major roads (primary/trunk/secondary by default)
    extracted from the real OSM graph.

Forbidden geometry (water, waterways, coastline) is fetched once via osmnx and
cached as a GeoPackage in `data/city_water/<city>.gpkg`.
s
Geometry units throughout the public API are metres (local UTM CRS chosen by
`osmnx.project_graph`). Defaults: `min_length=50` m, `min_angle=π/4`.
"""

from __future__ import annotations

import os
import sys
import warnings
from typing import Optional, Sequence

import geopandas as gpd
import numpy as np
import osmnx as ox
import pandas as pd
from scipy.spatial import cKDTree
from shapely.geometry import LineString, MultiLineString, Point, Polygon
from shapely.ops import unary_union
from shapely.prepared import prep
from shapely.strtree import STRtree

sys.path.insert(0, os.path.dirname(__file__))
from space_filling import (  # noqa: E402  reuse hot inner loops
    _connection_ok,
    _find_best,
    _segment_blocks,
    _min_incident_angle,
    EPS,
)


def nkey(p, decimals=6):
    """Coordinate key — 6 decimals (sub-mm) is robust to float64 roundtrip."""
    return (round(float(p[0]), decimals), round(float(p[1]), decimals))

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
GRAPH_DIR = os.path.join(REPO_ROOT, 'data', 'city_graphs')
WATER_DIR = os.path.join(REPO_ROOT, 'data', 'city_water')

WATER_TAGS = {
    'natural': ['water', 'coastline', 'bay', 'strait'],
    'waterway': True,
    'water': True,
    'landuse': ['reservoir', 'basin'],
}

DEFAULT_HIGHWAY_TAGS = ('motorway', 'trunk', 'primary', 'secondary')


# --------------------------------------------------------------------------- #
# Area + water fetch                                                          #
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


def _consolidate_water(gdf: gpd.GeoDataFrame, *, min_area_m2: float = 500.0,
                       simplify_tol_m: float = 5.0,
                       waterway_buffer_m: float = 3.0) -> gpd.GeoDataFrame:
    """Merge raw OSM water features into a small set of (multi)polygons.

    * Buffers waterway/coastline lines so they have non-zero width.
    * Drops polygons below `min_area_m2`.
    * Simplifies remaining polygons by `simplify_tol_m`.
    * Unions everything and returns a GeoDataFrame with one row per polygon.
    """
    if len(gdf) == 0:
        return gdf
    polys = []
    for geom in gdf.geometry:
        if geom is None or geom.is_empty:
            continue
        gt = geom.geom_type
        if gt in ('Polygon', 'MultiPolygon'):
            polys.append(geom)
        elif gt in ('LineString', 'MultiLineString'):
            polys.append(geom.buffer(waterway_buffer_m))
    if not polys:
        return gdf.iloc[0:0]
    merged = unary_union(polys)
    # split MultiPolygon into parts, filter small ones, simplify
    parts = list(merged.geoms) if merged.geom_type == 'MultiPolygon' else [merged]
    keep = [p.simplify(simplify_tol_m, preserve_topology=True)
            for p in parts if p.area >= min_area_m2]
    return gpd.GeoDataFrame({'geometry': keep}, geometry='geometry', crs=gdf.crs)


def fetch_water(area_poly_proj: Polygon, crs_proj, *, cache_path: Optional[str],
                min_area_m2: float = 500.0, simplify_tol_m: float = 5.0) -> gpd.GeoDataFrame:
    """Fetch water-like OSM features inside the area polygon.

    Returns a GeoDataFrame in `crs_proj` containing polygon and line geometries
    that mark forbidden boundaries. Cached to `cache_path` if given.
    """
    if cache_path and os.path.exists(cache_path):
        return gpd.read_file(cache_path)

    # osmnx wants WGS84 polygon for features_from_polygon
    area_wgs = gpd.GeoSeries([area_poly_proj], crs=crs_proj).to_crs(4326).iloc[0]
    try:
        feats = ox.features_from_polygon(area_wgs, tags=WATER_TAGS)
    except Exception as e:
        warnings.warn(f'water fetch failed: {e}; assuming no water')
        feats = gpd.GeoDataFrame(geometry=[], crs=4326)

    if len(feats) == 0:
        out = gpd.GeoDataFrame(geometry=[], crs=crs_proj)
    else:
        out = feats.to_crs(crs_proj)[['geometry']].copy()
        out = out[out.geometry.intersects(area_poly_proj)].copy()
        out['geometry'] = out.geometry.intersection(area_poly_proj)
        out = out[~out.geometry.is_empty]
        out = _consolidate_water(out, min_area_m2=min_area_m2,
                                 simplify_tol_m=simplify_tol_m)

    if cache_path:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        if len(out):
            out.to_file(cache_path, driver='GPKG')
        else:
            # write an empty marker so we don't refetch
            gpd.GeoDataFrame({'geometry': [area_poly_proj.centroid]},
                             geometry='geometry', crs=crs_proj).to_file(
                cache_path, driver='GPKG')
    return out


def load_or_fetch_water(city: str, G_proj, area_poly_proj: Polygon,
                       use_cache: bool = True) -> gpd.GeoDataFrame:
    cache = os.path.join(WATER_DIR, _city_safe_name(city) + '.gpkg') if use_cache else None
    return fetch_water(area_poly_proj, G_proj.graph['crs'], cache_path=cache)


# --------------------------------------------------------------------------- #
# Blocking segments (area boundary + water)                                   #
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


def _linestring_segments(line):
    coords = list(line.coords)
    for i in range(len(coords) - 1):
        yield (coords[i][:2], coords[i + 1][:2])


def boundary_segments(area_poly: Polygon) -> list:
    """Non-splittable blocking segments for the outer area boundary only.

    Water blocking is enforced separately via shapely (see `_make_water_check`)
    because OSM water can have tens of thousands of vertices and an O(N)
    crossing test against all of them is too slow.
    """
    return list(_polygon_boundary_segments(area_poly))


def _make_water_check(water_gdf: gpd.GeoDataFrame):
    """Return a `crosses_water(a, b) -> bool` predicate.

    Uses a prepared union of all water polygons + `line.difference(union)` so
    that a segment which has *any* portion strictly inside water is rejected,
    even when individual-polygon tests would only see touch points.
    Returns `None` if there is no water.
    """
    if water_gdf is None or len(water_gdf) == 0:
        return None
    geoms = [g for g in water_gdf.geometry if g is not None and not g.is_empty]
    if not geoms:
        return None
    water_union = unary_union(geoms)
    prep_union = prep(water_union)
    tol = 1e-3  # 1 mm — well below any meaningful street resolution

    def crosses_water(a, b) -> bool:
        ln = LineString([a, b])
        if not prep_union.intersects(ln):
            return False
        # Any part of the line clipped away by water means it crossed water.
        return ln.difference(water_union).length < ln.length - tol

    return crosses_water


# --------------------------------------------------------------------------- #
# Initialisation modes                                                        #
# --------------------------------------------------------------------------- #

def _land_polygon(area_poly: Polygon, water_gdf: gpd.GeoDataFrame):
    if water_gdf is None or len(water_gdf) == 0:
        return area_poly
    water_union = unary_union(list(water_gdf.geometry))
    land = area_poly.difference(water_union)
    return land if not land.is_empty else area_poly


def init_spokes(area_poly: Polygon, land_poly, n_spokes: int = 8,
                crosses_water=None) -> list:
    """Centroid → boundary spokes, clipped at the first water intersection."""
    c = land_poly.representative_point()
    cx, cy = c.x, c.y

    boundary = area_poly.exterior
    out = []
    for i in range(n_spokes):
        t = i / n_spokes
        p = boundary.interpolate(t, normalized=True)
        if crosses_water is None or not crosses_water((cx, cy), (p.x, p.y)):
            out.append(((cx, cy), (p.x, p.y)))
            continue
        # Clip the spoke at the first water boundary intersection: find the
        # furthest point along centroid→p that stays on land.
        ln = LineString([(cx, cy), (p.x, p.y)])
        clipped = ln.intersection(land_poly)
        if clipped.is_empty:
            continue
        if clipped.geom_type == 'MultiLineString':
            # take the piece containing the centroid (if any), else the longest
            pieces = list(clipped.geoms)
            near = [g for g in pieces if g.distance(c) < 1e-3]
            ln_use = near[0] if near else max(pieces, key=lambda g: g.length)
        else:
            ln_use = clipped
        if ln_use.length < 1.0:
            continue
        coords = list(ln_use.coords)
        out.append((coords[0][:2], coords[-1][:2]))
    return out


def init_from_roads(G_proj, area_poly: Polygon, land_poly,
                    highway_tags: Sequence[str] = DEFAULT_HIGHWAY_TAGS) -> list:
    """Edges from the projected OSM graph whose highway tag is in `highway_tags`,
    clipped to land (= area minus water)."""
    out = []
    for u, v, data in G_proj.edges(data=True):
        hw = data.get('highway')
        if isinstance(hw, list):
            hit = any(h in highway_tags for h in hw)
        else:
            hit = hw in highway_tags
        if not hit:
            continue
        geom = data.get('geometry')
        if geom is None:
            ux, uy = G_proj.nodes[u]['x'], G_proj.nodes[u]['y']
            vx, vy = G_proj.nodes[v]['x'], G_proj.nodes[v]['y']
            geom = LineString([(ux, uy), (vx, vy)])
        if not geom.intersects(land_poly):
            continue
        clipped = geom.intersection(land_poly)
        if clipped.geom_type == 'LineString':
            lns = [clipped]
        elif clipped.geom_type == 'MultiLineString':
            lns = list(clipped.geoms)
        else:
            continue
        for ln in lns:
            for a, b in _linestring_segments(ln):
                out.append((a, b))
    return out


# --------------------------------------------------------------------------- #
# Main simulation (real-area variant)                                         #
# --------------------------------------------------------------------------- #

def _find_best_water_safe(mx, my, candidates, idxs, seg_arr, model, min_angle,
                          check_water):
    for idx in idxs:
        tx, ty = candidates[idx]
        if not _connection_ok(mx, my, tx, ty, seg_arr, model, min_angle):
            continue
        if not check_water((float(tx), float(ty))):
            continue
        return (float(tx), float(ty))
    return None


def _run(seg_init, blocking_segs, area_poly, *, crosses_water=None,
         model=1, min_length=50.0, max_iter=200_000, seed=42,
         min_angle=np.pi / 4, area_coeff=2500.0, area_scale=None,
         kdtree_rebuild_every=256, progress_every=0, verbose=True):
    """Core loop. `seg_init` are splittable seeds; `blocking_segs` are not."""
    rng = np.random.default_rng(seed)

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
    # for the network centre (works whether seeds came from spokes or roads).
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
            return node_arr[gi[mask]]

        def _check(b):
            return crosses_water is None or not crosses_water((midx, midy), b)

        cand = _ordered(30)
        if len(cand) == 0:
            continue
        best = _find_best_water_safe(midx, midy, cand, range(len(cand)), sa,
                                     model, min_angle, _check)
        if best is None and tree_n > 30:  # widen search to all nodes
            cand_all = _ordered(tree_n)
            best = _find_best_water_safe(midx, midy, cand_all, range(len(cand_all)),
                                         sa, model, min_angle, _check)

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


def run_real_model(city: str, *, init: str = 'spokes', model: int = 1,
                   min_length: float = 50.0, min_angle: float = np.pi / 4,
                   max_iter: int = 200_000, seed: int = 42,
                   n_spokes: int = 8,
                   highway_tags: Sequence[str] = DEFAULT_HIGHWAY_TAGS,
                   area_coeff: float = 2500.0,
                   area_scale: Optional[float] = None,
                   kdtree_rebuild_every: int = 256,
                   use_water_cache: bool = True,
                   verbose: bool = True):
    """Run the real-area space-filling model for `city`.

    Returns
    -------
    dict with keys:
        segments      : list of ((x1,y1),(x2,y2)) in the local UTM CRS.
        area_polygon  : shapely Polygon (UTM metres).
        water         : GeoDataFrame of water features (UTM metres).
        crs           : pyproj CRS of the UTM projection.
    """
    if init not in ('spokes', 'roads'):
        raise ValueError(f'init must be "spokes" or "roads", got {init!r}')

    G_proj = load_city_graph(city)
    crs = G_proj.graph['crs']
    area = area_polygon_from_graph(G_proj)
    water = load_or_fetch_water(city, G_proj, area, use_cache=use_water_cache)
    land = _land_polygon(area, water)
    crosses_water = _make_water_check(water)

    if init == 'spokes':
        seeds = init_spokes(area, land, n_spokes=n_spokes,
                            crosses_water=crosses_water)
    else:
        seeds = init_from_roads(G_proj, area, land, highway_tags=highway_tags)
        if not seeds:
            warnings.warn('no major roads found in area; falling back to spokes')
            seeds = init_spokes(area, land, n_spokes=n_spokes,
                                crosses_water=crosses_water)

    blockers = boundary_segments(area)

    segments = _run(seeds, blockers, area,
                    crosses_water=crosses_water,
                    model=model, min_length=min_length, max_iter=max_iter,
                    seed=seed, min_angle=min_angle,
                    area_coeff=area_coeff, area_scale=area_scale,
                    kdtree_rebuild_every=kdtree_rebuild_every, verbose=verbose)

    return {
        'segments': segments,
        'area_polygon': area,
        'land_polygon': land,
        'water': water,
        'crs': crs,
        'init': init,
        'model': model,
    }


def run_circular_model(city: Optional[str] = None, *, area_m2: Optional[float] = None,
                       model: int = 1, min_length: float = 50.0,
                       min_angle: float = np.pi / 4, max_iter: int = 200_000,
                       seed: int = 42, n_spokes: int = 8,
                       area_coeff: float = 2500.0,
                       area_scale: Optional[float] = None,
                       kdtree_rebuild_every: int = 256,
                       progress_every: int = 1000,
                       crs=None, verbose: bool = True):
    """Run the space-filling model on a plain circular domain (no water).

    The circle's area matches the convex-hull area of `city` (if given) or the
    explicit `area_m2`. Useful for comparing model parameters in metres against
    real-city parameters without confounds from coastline geometry.
    """
    if area_m2 is None:
        if city is None:
            raise ValueError('Provide either `city` or `area_m2`')
        G_proj = load_city_graph(city)
        area_m2 = float(area_polygon_from_graph(G_proj).area)
        if crs is None:
            crs = G_proj.graph['crs']

    area = circular_area_polygon(area_m2)
    land = area
    seeds = init_spokes(area, land, n_spokes=n_spokes, crosses_water=None)
    blockers = boundary_segments(area)

    segments = _run(seeds, blockers, area, crosses_water=None,
                    model=model, min_length=min_length, max_iter=max_iter,
                    seed=seed, min_angle=min_angle,
                    area_coeff=area_coeff, area_scale=area_scale,
                    kdtree_rebuild_every=kdtree_rebuild_every,
                    progress_every=progress_every, verbose=verbose)

    return {
        'segments': segments,
        'area_polygon': area,
        'land_polygon': land,
        'water': gpd.GeoDataFrame(geometry=[], crs=crs),
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
    city = sys.argv[1] if len(sys.argv) > 1 else 'Copenhagen,_Denmark'
    print(f'=== Real space-filling model: {city} ===')
    for init in ('spokes', 'roads'):
        for model in (1, 2):
            print(f'\n[{init}, model {model}]')
            res = run_real_model(city, init=init, model=model,
                                 min_length=80.0, max_iter=20_000, seed=42)
            gdf = to_geodataframe(res, model_id=model)
            out = f'real_{_city_safe_name(city)}_{init}_m{model}.gpkg'
            gdf.to_file(out, driver='GPKG')
            print(f'  wrote {out} ({len(gdf)} segments)')
