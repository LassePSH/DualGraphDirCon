#!/usr/bin/env python3
"""
Space-filling street network model.

Recreates the two models described in the urban street pattern paper:

Model 1 (basic): At each step, pick a segment proportional to length, split at
  midpoint, connect midpoint to the closest visible intersection.

Model 2 (biased): Same as Model 1 but enforces:
  - Polygon area constraint: A(r) > 0.05 * exp(-1/r), larger blocks near periphery
  - Angle constraint: all angles at intersections must be > pi/4

Both start from a unit square divided by 4 lines (horizontal, vertical, two
diagonals) all crossing at the centre.

Requires: numpy, shapely, geopandas
"""

import numpy as np
from shapely.geometry import LineString
from collections import defaultdict
from scipy.spatial import cKDTree
import geopandas as gpd

# --------------------------------------------------------------------------- #
# Parameters                                                                    #
# --------------------------------------------------------------------------- #

MIN_LENGTH = 0.05       # Minimum segment length; stop when nothing is longer
CENTER     = (0.5, 0.5) # Centre of the unit square
EPS        = 1e-9       # Floating-point tolerance

# --------------------------------------------------------------------------- #
# Geometry helpers                                                            #
# --------------------------------------------------------------------------- #

def dist(p1, p2):
    return np.hypot(p1[0] - p2[0], p1[1] - p2[1])


def midpt(p1, p2):
    return ((p1[0] + p2[0]) / 2.0, (p1[1] + p2[1]) / 2.0)


def seg_len(seg):
    return dist(seg[0], seg[1])


def properly_crosses(p1, p2, p3, p4):
    """Return True if segment p1-p2 crosses segment p3-p4 at an interior point."""
    s1 = LineString([p1, p2])
    s2 = LineString([p3, p4])
    inter = s1.intersection(s2)
    if inter.is_empty:
        return False
    if inter.geom_type == 'Point':
        x, y = inter.x, inter.y
        for ep in (p1, p2, p3, p4):
            if abs(x - ep[0]) < EPS and abs(y - ep[1]) < EPS:
                return False   # touching at an endpoint only
        return True
    return True  # line overlap


def is_visible(p, q, segments):
    """Can p see q without any segment (not incident to p or q) blocking the view?"""
    for a, b in segments:
        if dist(a, p) < EPS or dist(b, p) < EPS:
            continue
        if dist(a, q) < EPS or dist(b, q) < EPS:
            continue
        if properly_crosses(p, q, a, b):
            return False
    return True


def min_angle_at(node, new_vec, segments):
    """Minimum angle (radians) between new_vec and every existing edge at node."""
    nv = np.asarray(new_vec, dtype=float)
    nv_norm = np.linalg.norm(nv)
    min_ang = np.pi
    for a, b in segments:
        if dist(a, node) < EPS:
            other = b
        elif dist(b, node) < EPS:
            other = a
        else:
            continue
        ev = np.asarray(other, dtype=float) - np.asarray(node, dtype=float)
        ev_norm = np.linalg.norm(ev)
        if ev_norm < EPS or nv_norm < EPS:
            continue
        cos_a = np.clip(np.dot(nv, ev) / (nv_norm * ev_norm), -1.0, 1.0)
        min_ang = min(min_ang, np.arccos(cos_a))
    return min_ang


def angles_ok(m, target, segments, min_angle=np.pi / 4):
    """Check angle constraint at both endpoints of the proposed new edge m->target."""
    vec = (target[0] - m[0], target[1] - m[1])
    if min_angle_at(m,      vec,                   segments) < min_angle:
        return False
    if min_angle_at(target, (-vec[0], -vec[1]),    segments) < min_angle:
        return False
    return True


# --------------------------------------------------------------------------- #
# Vectorized geometry helpers (used by run_model for speed)                    #
# --------------------------------------------------------------------------- #

def _is_visible_vec(p, q, seg_arr):
    """Vectorized visibility check against a numpy segment array (N, 2, 2)."""
    if len(seg_arr) == 0:
        return True
    px, py = p
    qx, qy = q
    ax = seg_arr[:, 0, 0];  ay = seg_arr[:, 0, 1]
    bx = seg_arr[:, 1, 0];  by = seg_arr[:, 1, 1]

    # Exclude segments incident to p or q
    inc = (((np.abs(ax - px) < EPS) & (np.abs(ay - py) < EPS)) |
           ((np.abs(bx - px) < EPS) & (np.abs(by - py) < EPS)) |
           ((np.abs(ax - qx) < EPS) & (np.abs(ay - qy) < EPS)) |
           ((np.abs(bx - qx) < EPS) & (np.abs(by - qy) < EPS)))
    ax, ay, bx, by = ax[~inc], ay[~inc], bx[~inc], by[~inc]
    if len(ax) == 0:
        return True

    # Proper crossing test via cross products
    d1 = (qx - px) * (ay - py) - (qy - py) * (ax - px)
    d2 = (qx - px) * (by - py) - (qy - py) * (bx - px)
    d3 = (bx - ax) * (py - ay) - (by - ay) * (px - ax)
    d4 = (bx - ax) * (qy - ay) - (by - ay) * (qx - ax)
    return not np.any((d1 * d2 < 0) & (d3 * d4 < 0))


def _min_angle_at_vec(node, new_vec, seg_arr):
    """Vectorized minimum angle at node against a numpy segment array (N, 2, 2)."""
    if len(seg_arr) == 0:
        return np.pi
    nx, ny = node
    ax = seg_arr[:, 0, 0];  ay = seg_arr[:, 0, 1]
    bx = seg_arr[:, 1, 0];  by = seg_arr[:, 1, 1]

    inc_a = (np.abs(ax - nx) < EPS) & (np.abs(ay - ny) < EPS)
    inc_b = (np.abs(bx - nx) < EPS) & (np.abs(by - ny) < EPS)
    incident = inc_a | inc_b
    if not incident.any():
        return np.pi

    ox = np.where(inc_a[incident], bx[incident], ax[incident])
    oy = np.where(inc_a[incident], by[incident], ay[incident])
    evx, evy = ox - nx, oy - ny
    ev_norms = np.hypot(evx, evy)

    nv = np.asarray(new_vec, dtype=float)
    nv_norm = np.linalg.norm(nv)
    if nv_norm < EPS:
        return np.pi

    valid = ev_norms > EPS
    if not valid.any():
        return np.pi

    cos_a = np.clip(
        (nv[0] * evx[valid] + nv[1] * evy[valid]) / (nv_norm * ev_norms[valid]),
        -1.0, 1.0,
    )
    return float(np.arccos(cos_a).min())


def _angles_ok_vec(m, target, seg_arr, min_angle=np.pi / 4):
    vec = (target[0] - m[0], target[1] - m[1])
    if _min_angle_at_vec(m,      vec,                  seg_arr) < min_angle:
        return False
    if _min_angle_at_vec(target, (-vec[0], -vec[1]),   seg_arr) < min_angle:
        return False
    return True


def area_ok(seg, mid, area_coeff=0.05, area_scale=1.0):
    """
    Model 2 area constraint: the approximate polygon area after splitting
    must satisfy A(r) > area_coeff * exp(-area_scale/r), where r is distance
    from centre. A polygon's area is approximated as (half-segment-length)^2.
    """
    r = dist(mid, CENTER)
    if r < EPS:
        return True
    min_area = area_coeff * np.exp(-area_scale / r)
    approx_area = (seg_len(seg) / 2.0) ** 2
    return approx_area >= min_area


# --------------------------------------------------------------------------- #
# Node bookkeeping                                                               #
# --------------------------------------------------------------------------- #

def nkey(p, decimals=8):
    """Hashable key for a point, rounded to avoid floating-point duplicates."""
    return (round(p[0], decimals), round(p[1], decimals))


def unique_nodes(segments):
    """Return a list of all unique nodes appearing in segments."""
    seen = {}
    for a, b in segments:
        seen[nkey(a)] = a
        seen[nkey(b)] = b
    return list(seen.values())


# --------------------------------------------------------------------------- #
# Initialisation                                                                 #
# --------------------------------------------------------------------------- #

def initialize():
    """
    Unit square divided by 4 lines (horizontal, vertical, two diagonals)
    all passing through the centre -> 8 initial segments.
    """
    c = CENTER
    line_ends = [
        ((0.0, 0.5), (1.0, 0.5)),   # horizontal
        ((0.5, 0.0), (0.5, 1.0)),   # vertical
        ((0.0, 0.0), (1.0, 1.0)),   # diagonal
        ((1.0, 0.0), (0.0, 1.0)),   # anti-diagonal
    ]
    segs = []
    for p1, p2 in line_ends:
        segs.append((p1, c))
        segs.append((c, p2))
    return segs


# --------------------------------------------------------------------------- #
# Main simulation                                                                #
# --------------------------------------------------------------------------- #

def run_model(model=1, min_length=MIN_LENGTH, max_iter=200_000, seed=42,
              min_angle=np.pi / 4, area_coeff=0.05, area_scale=1.0):
    """
    Run the space-filling model.

    Parameters
    ----------
    model      : 1 (basic) or 2 (with area + angle biases)
    min_length : segments shorter than this are never split
    max_iter   : hard iteration cap
    seed       : random seed
    min_angle  : (model 2) minimum allowed angle at intersections in radians
    area_coeff : (model 2) coefficient in area constraint: area_coeff * exp(-area_scale/r)
    area_scale : (model 2) scale factor in area constraint exponent
    """
    np.random.seed(seed)
    init = initialize()

    # Pre-allocate segment array with capacity doubling (avoids O(N) Python→C
    # conversions on every iteration).
    capacity = max(len(init) * 4, 64)
    seg_arr = np.empty((capacity, 2, 2), dtype=np.float64)
    n = len(init)
    for i, (a, b) in enumerate(init):
        seg_arr[i, 0] = a
        seg_arr[i, 1] = b

    # Incremental node tracking: nkey → index in node_arr
    node_keys  = {}   # nkey(pt) -> row index in node_arr
    node_cap   = capacity * 2
    node_arr   = np.empty((node_cap, 2), dtype=np.float64)
    n_nodes    = 0
    for a, b in init:
        for pt in (a, b):
            k = nkey(pt)
            if k not in node_keys:
                if n_nodes >= node_cap:
                    node_cap *= 2
                    node_arr = np.resize(node_arr, (node_cap, 2))
                node_arr[n_nodes] = pt
                node_keys[k] = n_nodes
                n_nodes += 1

    def _add_node(pt):
        nonlocal n_nodes, node_cap, node_arr
        k = nkey(pt)
        if k not in node_keys:
            if n_nodes >= node_cap:
                node_cap *= 2
                node_arr = np.resize(node_arr, (node_cap, 2))
            node_arr[n_nodes] = pt
            node_keys[k] = n_nodes
            n_nodes += 1
        return node_keys[k]

    def _ensure_seg_cap(needed):
        nonlocal capacity, seg_arr
        if needed > capacity:
            capacity = max(capacity * 2, needed)
            new = np.empty((capacity, 2, 2), dtype=np.float64)
            new[:n] = seg_arr[:n]
            seg_arr = new

    for it in range(max_iter):
        sa = seg_arr[:n]   # live view — no copy

        lengths = np.hypot(
            sa[:, 1, 0] - sa[:, 0, 0],
            sa[:, 1, 1] - sa[:, 0, 1],
        )

        eligible_mask = lengths > min_length

        # Model 2: vectorised area constraint
        if model == 2:
            mids = (sa[:, 0] + sa[:, 1]) * 0.5
            r = np.hypot(mids[:, 0] - CENTER[0], mids[:, 1] - CENTER[1])
            with np.errstate(divide='ignore', invalid='ignore'):
                min_area_arr = np.where(
                    r > EPS, area_coeff * np.exp(-area_scale / r), 0.0
                )
            eligible_mask &= (lengths * 0.5) ** 2 >= min_area_arr

        eligible_indices = np.where(eligible_mask)[0]
        if len(eligible_indices) == 0:
            print(f"  Model {model}: converged at iteration {it} with {n} segments")
            break

        # Sample proportional to length
        elig_lengths = lengths[eligible_indices]
        probs = elig_lengths / elig_lengths.sum()
        choice_idx = int(np.random.choice(len(eligible_indices), p=probs))
        seg_i = int(eligible_indices[choice_idx])

        a = tuple(sa[seg_i, 0])
        b = tuple(sa[seg_i, 1])
        m = midpt(a, b)

        # Swap-remove seg_i: O(1), no shifting
        n -= 1
        seg_arr[seg_i] = seg_arr[n]

        # Append (a,m) and (m,b)
        _ensure_seg_cap(n + 2)
        seg_arr[n,   0] = a;  seg_arr[n,   1] = m;  n += 1
        seg_arr[n,   0] = m;  seg_arr[n,   1] = b;  n += 1

        # Register new node m (a and b already tracked)
        m_idx = _add_node(m)
        a_idx = node_keys[nkey(a)]
        b_idx = node_keys[nkey(b)]

        # Build candidate array: all nodes except m, a, b
        all_nodes = node_arr[:n_nodes]
        excl = np.array([m_idx, a_idx, b_idx])
        mask = np.ones(n_nodes, dtype=bool)
        mask[excl] = False
        candidates = all_nodes[mask]

        if len(candidates) == 0:
            continue

        sa = seg_arr[:n]   # refresh view after modifications

        # Query only the k nearest neighbours — avoids full sort (huge speedup).
        # Fall back to all candidates only if the initial batch yields nothing.
        K_INIT = min(30, len(candidates))
        tree = cKDTree(candidates)
        dists, indices = tree.query(m, k=K_INIT, workers=1)
        if K_INIT == 1:
            dists, indices = [float(dists)], [int(indices)]

        best_node = None
        for idx in indices:
            node = tuple(candidates[idx])
            if not _is_visible_vec(m, node, sa):
                continue
            if model == 2 and not _angles_ok_vec(m, node, sa, min_angle):
                continue
            best_node = node
            break

        # Fallback: expand to all candidates (rare)
        if best_node is None and K_INIT < len(candidates):
            _, indices_all = tree.query(m, k=len(candidates), workers=1)
            for idx in indices_all[K_INIT:]:
                node = tuple(candidates[idx])
                if not _is_visible_vec(m, node, sa):
                    continue
                if model == 2 and not _angles_ok_vec(m, node, sa, min_angle):
                    continue
                best_node = node
                break

        if best_node is not None:
            _ensure_seg_cap(n + 1)
            seg_arr[n, 0] = m
            seg_arr[n, 1] = best_node
            n += 1

    else:
        print(f"  Model {model}: hit max_iter={max_iter}, {n} segments")

    return [(tuple(seg_arr[i, 0]), tuple(seg_arr[i, 1])) for i in range(n)]


# --------------------------------------------------------------------------- #
# Shapefile export                                                               #
# --------------------------------------------------------------------------- #

def compute_degrees(segments):
    """Return a dict mapping each node key to its degree."""
    deg = defaultdict(int)
    for a, b in segments:
        deg[nkey(a)] += 1
        deg[nkey(b)] += 1
    return dict(deg)


def to_geodataframe(segments, model_id):
    """
    Convert a list of (point, point) segments to a GeoDataFrame.

    Attributes per segment:
      model     – model number (1 or 2)
      length    – Euclidean length of the segment
      deg_start – degree of the start node
      deg_end   – degree of the end node
    """
    degrees = compute_degrees(segments)
    records = []
    for a, b in segments:
        records.append({
            'geometry':  LineString([a, b]),
            'model':     model_id,
            'length':    seg_len((a, b)),
            'deg_start': degrees.get(nkey(a), 0),
            'deg_end':   degrees.get(nkey(b), 0),
        })
    return gpd.GeoDataFrame(records, crs='EPSG:4326')


def save_shapefile(segs, model_id, path):
    gdf = to_geodataframe(segs, model_id)
    gdf.to_file(path)
    print(f"  Saved {len(gdf)} segments → {path}")
    return gdf


def get_geodataframe(model=1, min_length=MIN_LENGTH, max_iter=200_000, seed=42,
                     min_angle=np.pi / 4, area_coeff=0.05, area_scale=1.0):
    """
    Run the model and return the result as a GeoDataFrame (no file I/O).

    Parameters
    ----------
    model      : 1 (basic) or 2 (with area + angle biases)
    min_length : segments shorter than this are never split
    max_iter   : hard iteration cap
    seed       : random seed
    min_angle  : (model 2) minimum allowed angle at intersections in radians
    area_coeff : (model 2) coefficient in area constraint: area_coeff * exp(-area_scale/r)
    area_scale : (model 2) scale factor in area constraint exponent
    """
    segs = run_model(model=model, min_length=min_length, max_iter=max_iter, seed=seed,
                     min_angle=min_angle, area_coeff=area_coeff, area_scale=area_scale)
    return to_geodataframe(segs, model_id=model)


# --------------------------------------------------------------------------- #
# Entry point                                                                    #
# --------------------------------------------------------------------------- #

if __name__ == '__main__':
    print("=== Space-Filling Street Network Model ===\n")

    print("Running Model 1 (basic)...")
    segs1 = run_model(model=1, seed=42)
    save_shapefile(segs1, model_id=1, path='model1.shp')

    print("Running Model 2 (with area & angle biases)...")
    segs2 = run_model(model=2, seed=42)
    save_shapefile(segs2, model_id=2, path='model2.shp')
