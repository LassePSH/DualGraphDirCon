import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from multiprocessing import Pool
from dgdc import get_dual_dir_con
from download_cities import load_city

CITY = 'Copenhagen, Denmark'
BASE = '/home/lpsha/s154446/fractality/data/city_degrees'

T_S = [5, 10, 50, 100]
A_S = [5, 10, 20, 30, 40]
SIMPLIFY = [True, False]


def folder_for(t_buf, a_thresh, s):
    suffix = '' if s else '_no_simplify'
    return f'{BASE}/t{int(t_buf)}_a{int(a_thresh)}{suffix}/'


def run_one(args):
    t_buf, a_thresh, s = args
    folder = folder_for(t_buf, a_thresh, s)
    out_file = folder + CITY + '.npz'
    if os.path.exists(out_file):
        return f't={t_buf} a={a_thresh} s={int(s)}: skipped (exists)'
    try:
        G = load_city(CITY)
        gdf_merged, _, _, _ = get_dual_dir_con(
            t_buffer=t_buf,
            a_threshold=a_thresh,
            data=G,
            simplify_roundabout=s,
            enforce_degree2=False,
        )
        os.makedirs(folder, exist_ok=True)
        np.savez(out_file, degree=np.array(gdf_merged.degree), length=np.array(gdf_merged.length))
        return f't={t_buf} a={a_thresh} s={int(s)}: success'
    except Exception as e:
        return f't={t_buf} a={a_thresh} s={int(s)}: failed: {e}'


if __name__ == '__main__':
    combos = [(t, a, s) for t in T_S for a in A_S for s in SIMPLIFY]
    print(f'Recomputing {CITY} across {len(combos)} combos')
    with Pool(20) as p:
        for msg in p.imap_unordered(run_one, combos):
            print(msg, flush=True)
    print('Done.')
