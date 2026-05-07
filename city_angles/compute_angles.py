"""
Compute the angle between every pair of intersecting streets for each cached
city graph in ``data/city_graphs/``.

For each city we build the dual graph (streets become nodes, shared
intersections become edges) using the same pipeline as ``dgdc.get_dual_dir_con``
up to the ``new_angles`` step, then dump the per-edge angle (degrees, acute
0-90) to ``city_angles/<city>.out`` as a 1-D ``np.savetxt`` array.
"""

import os
import sys
import argparse
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import osmnx as ox
import momepy
import geopandas as gpd

from dgdc.dual_conti import clean_chains, new_angles

GRAPH_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'city_graphs')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'city_angles')


def compute_angles(city_file, t_buffer):
    G = ox.load_graphml(os.path.join(GRAPH_DIR, city_file))

    shape_df = ox.graph_to_gdfs(ox.convert.to_undirected(G), nodes=False)
    shape_df = shape_df.set_crs("epsg:4326")
    shape_df = shape_df.to_crs(shape_df.estimate_utm_crs())

    if 'edgeUID' not in shape_df.columns:
        shape_df['edgeUID'] = shape_df.index

    shape_df = shape_df.reset_index().explode('geometry')
    u = shape_df.union_all()
    out = gpd.GeoDataFrame(
        geometry=gpd.GeoSeries(u, crs=shape_df.crs).explode()
    ).reset_index(drop=True)
    shape_exploded_df = out.sjoin(
        shape_df[['osmid', 'geometry', 'edgeUID']], how="left", predicate="intersects"
    )
    shape_exploded_df = shape_exploded_df.drop_duplicates(subset=['geometry']).reset_index(drop=True)
    shape_exploded_df['id'] = shape_exploded_df.index

    G_primal = momepy.gdf_to_nx(shape_exploded_df, approach="primal")
    G_primal = clean_chains(G_primal)
    _, lines = momepy.nx_to_gdf(G_primal)

    G_dual = momepy.gdf_to_nx(lines, approach='dual', multigraph=False, angles=False)
    G_dual = new_angles(G_dual, touch_buffer=t_buffer)

    return np.array([d['new_angle'] for _, _, d in G_dual.edges(data=True)])


def parallel(city_file):
    name = os.path.splitext(city_file)[0]
    out_path = os.path.join(OUT_DIR, name + '.out')
    if os.path.exists(out_path):
        print(name + ' skipped (exists)')
        return
    try:
        angles = compute_angles(city_file, t_buffer)
        np.savetxt(out_path, angles, delimiter=',')
        print(name + ' success')
    except Exception as e:
        print(f'{name} failed: {e}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute per-intersection street angles for each cached city graph.'
    )
    parser.add_argument('--t_buffer', type=float, default=50)
    parser.add_argument('--processes', type=int, default=20)
    args = parser.parse_args()

    t_buffer = args.t_buffer

    os.makedirs(OUT_DIR, exist_ok=True)

    cities = sorted(f for f in os.listdir(GRAPH_DIR) if f.endswith('.graphml'))
    print(f'Processing {len(cities)} cities with t_buffer={t_buffer} ...')

    with Pool(args.processes) as p:
        p.map(parallel, cities)
