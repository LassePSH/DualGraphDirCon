"""Compute degree distribution of the dual graph induced by COINS strokes.

Each COINS stroke becomes a dual-graph node; two strokes are connected if they
share an endpoint in the primal network. Degree of a stroke is the number of
distinct other strokes meeting it at either endpoint.

COINS angle_threshold is the minimum *interior* angle for two segments to be
grouped (180 = perfectly straight), so to mirror DGDC's max-deflection
a_threshold we pass 180 - a_threshold.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import argparse
import warnings
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
import numpy as np
import osmnx as ox
import momepy

from cities.download_cities import load_city


def compute_coins_dual_degrees(G, a_threshold=20):
    """Return (degree, length) numpy arrays, one entry per COINS stroke.

    Dual-graph construction: at every primal intersection node, all distinct
    stroke IDs in incident edges are mutually connected. Degree of a stroke is
    the number of distinct other strokes it meets at any shared intersection.

    Parameters
    ----------
    G : networkx graph from osmnx
    a_threshold : float
        Max deflection in degrees (DGDC convention). Translated to
        COINS angle_threshold = 180 - a_threshold.
    """
    G_proj = ox.project_graph(ox.convert.to_undirected(G))
    _, edges = ox.graph_to_gdfs(G_proj)
    edges = edges[['geometry']].reset_index(drop=True)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        coins = momepy.COINS(edges, angle_threshold=180 - a_threshold)
        stroke_attr = coins.stroke_attribute()
        strokes = coins.stroke_gdf()

    edges = edges.copy()
    edges['stroke'] = stroke_attr.values

    point_to_strokes = defaultdict(set)
    for sid, geom in zip(edges['stroke'], edges.geometry):
        coords = list(geom.coords)
        point_to_strokes[coords[0]].add(sid)
        point_to_strokes[coords[-1]].add(sid)

    neighbors = defaultdict(set)
    for sid_set in point_to_strokes.values():
        if len(sid_set) < 2:
            continue
        for sid in sid_set:
            neighbors[sid].update(sid_set)
    for sid in neighbors:
        neighbors[sid].discard(sid)

    sids = list(strokes.index)
    degree = np.array([len(neighbors[sid]) for sid in sids], dtype=int)
    length = strokes.geometry.length.to_numpy()
    return degree, length


def run_city(city, a_threshold, out_dir, skip_existing=True):
    safe = city.replace('/', '_')
    out_path = os.path.join(out_dir, safe + '.npz')
    if skip_existing and os.path.exists(out_path):
        print(f'{city}: already done, skipping')
        return
    try:
        G = load_city(city)
        degree, length = compute_coins_dual_degrees(G, a_threshold=a_threshold)
        np.savez(out_path, degree=degree, length=length)
        print(f'{city}: {len(degree)} strokes -> {out_path}')
    except Exception as e:
        print(f'{city}: FAILED ({type(e).__name__}: {e})')


def _worker(city, a_threshold, out_dir):
    run_city(city, a_threshold, out_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='COINS dual-graph degree distribution per city.')
    parser.add_argument('--city', type=str, default=None,
                        help='Single city (must be cached in data/city_graphs/). '
                             'If omitted, runs all cities in cities/cities.txt.')
    parser.add_argument('--a_threshold', type=float, default=20,
                        help='Max deflection in degrees (DGDC convention).')
    parser.add_argument('--processes', type=int, default=20)
    parser.add_argument('--out_dir', type=str, default=None)
    args = parser.parse_args()

    out_dir = args.out_dir or os.path.join(
        os.path.dirname(__file__), '..', 'data',
        f'city_degrees_coins/a{int(args.a_threshold)}'
    )
    os.makedirs(out_dir, exist_ok=True)

    if args.city is not None:
        run_city(args.city, args.a_threshold, out_dir)
    else:
        cities_file = os.path.join(os.path.dirname(__file__), '..', 'cities', 'cities.txt')
        with open(cities_file) as f:
            cities = [line.strip() for line in f if line.strip()]
        print(f'Running COINS on {len(cities)} cities with {args.processes} processes')
        worker = partial(_worker, a_threshold=args.a_threshold, out_dir=out_dir)
        with Pool(args.processes) as p:
            p.map(worker, cities)
