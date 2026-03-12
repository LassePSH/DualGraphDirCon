import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dgdc import get_dual_dir_con
import numpy as np
from multiprocessing import Pool
import osmnx as ox
import argparse
import os


def get_degree_seq(city, t_buffer, a_threshold):
    G = ox.graph.graph_from_place(city)
    gdf_merged, _, _, _ = get_dual_dir_con(t_buffer=t_buffer, a_threshold=a_threshold, data=G, enforce_degree2=False)
    degree_sequence = np.array(gdf_merged.degree)
    return degree_sequence

def parallel(city):
    try:
        degree_sequence = get_degree_seq(city, t_buffer, a_threshold)
        np.savetxt(path+city+'.out', degree_sequence, delimiter=','); del degree_sequence
        print(city + ' succes!')
    except:
        print(city + ' failed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute dual graph degree sequences for a set of cities using OSM street networks.',
        epilog='Example: python regime.py --t_buffer 10 --a_threshold 20'
    )
    parser.add_argument('--t_buffer', type=float, default=10)
    parser.add_argument('--a_threshold', type=float, default=20)
    args = parser.parse_args()

    t_buffer = args.t_buffer
    a_threshold = args.a_threshold

    cities_file = os.path.join(os.path.dirname(__file__), 'cities.txt')
    with open(cities_file) as f:
        cities = [line.strip() for line in f if line.strip()]

    print("Number of cities:", len(cities))

    path = f'/home/lpsha/s154446/fractality/data/city_degrees/t{int(t_buffer)}_a{int(a_threshold)}/'

    with Pool(20) as p:
        p.map(parallel, cities)