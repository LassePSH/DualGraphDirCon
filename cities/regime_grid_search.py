import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dgdc import get_dual_dir_con
import numpy as np
from multiprocessing import Pool
from download_cities import load_city


t_buffer = None
a_threshold = None
path = None


def get_degree_length(city, t_buf, a_thresh,s):
    G = load_city(city)
    gdf_merged, _, _, _ = get_dual_dir_con(t_buffer=t_buf, 
                                           a_threshold=a_thresh, 
                                           data=G, 
                                           simplify_roundabout=s, 
                                           enforce_degree2=False)
    
    return np.array(gdf_merged.degree), np.array(gdf_merged.length)

def parallel(city):
    out_file = path + city + '.npz'
    if os.path.exists(out_file):
        print(city + ' skipped (exists)')
        return
    try:
        degree, length = get_degree_length(city, t_buffer, a_threshold, s)
        np.savez(out_file, degree=degree, length=length)
        print(city + ' success!')
    except Exception as e:
        print(city + f' failed: {e}')


if __name__ == '__main__':
    t_s = [5, 10, 50, 100]
    a_s = [5, 10, 20, 30, 40]
    simplify = [True, False]

    cities_file = os.path.join(os.path.dirname(__file__), 'cities.txt')
    with open(cities_file) as f:
        cities = [line.strip() for line in f if line.strip()]

    print(f"Number of cities: {len(cities)}")
    print(f"Grid: t_buffer={t_s}, a_threshold={a_s}  ({len(t_s)*len(a_s)} combinations)\n")

    for t_buf in t_s:
        for a_thresh in a_s:
            for s in simplify:
                t_buffer = t_buf
                a_threshold = a_thresh

                if s:
                    path = f'/home/lpsha/s154446/fractality/data/city_degrees/t{int(t_buf)}_a{int(a_thresh)}/'
                else: 
                    path = f'/home/lpsha/s154446/fractality/data/city_degrees/t{int(t_buf)}_a{int(a_thresh)}_no_simplify/'

                if os.path.exists(path):
                    print(f't={t_buf} a={a_thresh}: exist, skipping folder')
                    continue

                os.makedirs(path, exist_ok=True)
                print(f't={t_buf} a={a_thresh} simplify={int(s)}:  running...')

                with Pool(50) as p:
                    p.map(parallel, cities)

