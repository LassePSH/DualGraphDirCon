"""
Compute per-edge compass bearings for each cached city graph, using two
methods following Boeing (2019):

Method 1 (simplified topology): nodes only at intersections/dead-ends, edges
    span whole street segments. Bearing of edge u->v is the compass heading
    from u to v, and its reciprocal (bearing + 180) mod 360 is also recorded,
    since a street centerline points in both directions.

Method 2 (unsimplified topology): edges are the raw straight-line pieces
    making up each simplified edge's ``LineString`` geometry (these are OSM's
    original inter-node segments preserved by ``osmnx`` during simplification).
    Each piece's bearing is emitted with its own length as a weight to correct
    for very short pieces that approximate curving streets.

Self-loops have undefined bearing and are skipped in both methods.

Both methods read only from ``data/city_graphs/`` (populated by
``cities/download_cities.py``). No network access.

Outputs one ``.npz`` per city to ``data/city_orientations/<city>.npz`` with:
    method1_bearings : 1-D float array, degrees in [0, 360)
    method2_bearings : 1-D float array, degrees in [0, 360)
    method2_weights  : 1-D float array, segment lengths (metres) aligned with
                       ``method2_bearings``
"""

import os
import sys
import argparse
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import osmnx as ox

GRAPH_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'city_graphs')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'city_orientations')


def _city_name_from_file(city_file):
    return os.path.splitext(city_file)[0]


def _method1_bearings(G):
    """Endpoint-to-endpoint bearing per simplified edge, plus reciprocal."""
    G = ox.convert.to_undirected(G)
    G = ox.bearing.add_edge_bearings(G)
    bearings = []
    for u, v, data in G.edges(data=True):
        if u == v:
            continue
        b = data.get('bearing')
        if b is None or np.isnan(b):
            continue
        bearings.append(b)
        bearings.append((b + 180.0) % 360.0)
    return np.asarray(bearings, dtype=float)


def _compass_bearing(lat1, lon1, lat2, lon2):
    """Great-circle initial bearing from (lat1,lon1) to (lat2,lon2) in degrees."""
    phi1 = np.radians(lat1)
    phi2 = np.radians(lat2)
    dlam = np.radians(lon2 - lon1)
    y = np.sin(dlam) * np.cos(phi2)
    x = np.cos(phi1) * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.cos(dlam)
    return (np.degrees(np.arctan2(y, x)) + 360.0) % 360.0


def _method2_bearings(G):
    """Bearing per raw coordinate pair inside each edge's LineString geometry,
    weighted by the metric length of that piece."""
    G_u = ox.convert.to_undirected(G)
    bearings, lengths = [], []
    for u, v, data in G_u.edges(data=True):
        if u == v:
            continue
        geom = data.get('geometry')
        if geom is None:
            continue
        coords = np.asarray(geom.coords)
        if len(coords) < 2:
            continue
        lons, lats = coords[:, 0], coords[:, 1]
        b = _compass_bearing(lats[:-1], lons[:-1], lats[1:], lons[1:])

        # metric length per piece via equirectangular approx at local latitude
        R = 6371000.0
        phi_m = np.radians(0.5 * (lats[:-1] + lats[1:]))
        dlat = np.radians(lats[1:] - lats[:-1])
        dlon = np.radians(lons[1:] - lons[:-1])
        L = R * np.hypot(dlat, dlon * np.cos(phi_m))

        mask = np.isfinite(b) & (L > 0)
        b, L = b[mask], L[mask]

        bearings.append(b)
        lengths.append(L)
        bearings.append((b + 180.0) % 360.0)
        lengths.append(L)

    if not bearings:
        return np.empty(0), np.empty(0)
    return np.concatenate(bearings), np.concatenate(lengths)


def compute_orientations(city_file):
    G = ox.load_graphml(os.path.join(GRAPH_DIR, city_file))
    method1 = _method1_bearings(G)
    method2, weights = _method2_bearings(G)
    return method1, method2, weights


def parallel(city_file):
    name = _city_name_from_file(city_file)
    out_path = os.path.join(OUT_DIR, name + '.npz')
    if os.path.exists(out_path):
        print(name + ' skipped (exists)')
        return
    try:
        m1, m2, w = compute_orientations(city_file)
        np.savez(out_path, method1_bearings=m1, method2_bearings=m2, method2_weights=w)
        print(f'{name} success (m1={len(m1)}, m2={len(m2)})')
    except Exception as e:
        print(f'{name} failed: {e}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute per-edge compass bearings for each cached city graph.'
    )
    parser.add_argument('--processes', type=int, default=20)
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    cities = sorted(f for f in os.listdir(GRAPH_DIR) if f.endswith('.graphml'))
    print(f'Processing {len(cities)} cities ...')

    with Pool(args.processes) as p:
        p.map(parallel, cities)
