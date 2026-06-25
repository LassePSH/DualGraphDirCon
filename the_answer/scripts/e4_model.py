"""E4: Space-filling model decomposition, with and without noise decoration.

Runs the growth model once (per model id / seed), then decorates the SAME
backbone with add_noise_streets at several noise_frac levels. Each network is
pushed through the identical instrumented DGDC decomposition used for cities,
plus per-street noise tags (a street is 'noise' if its segments are not part
of the pre-noise backbone).
"""
import os
import sys
import time
import warnings
warnings.filterwarnings('ignore')

import numpy as np

sys.path.insert(0, '/home/lpsha/s154446/fractality')
sys.path.insert(0, '/home/lpsha/s154446/fractality/space_filling_model')
sys.path.insert(0, '/home/lpsha/s154446/fractality/the_answer/scripts')

OUT_DIR = '/home/lpsha/s154446/fractality/the_answer/results/model_decomp'
MIN_LENGTH = 0.004
A_THRESHOLD = 20
T_BUFFER = 0.0005  # irrelevant for straight segments; must be < min_length


def segments_to_gdf(segs):
    import geopandas as gpd
    from shapely.geometry import LineString
    return gpd.GeoDataFrame(geometry=[LineString([a, b]) for a, b in segs])


def decompose_segments(segs, base_union=None):
    """Run the city decomposition on a model segment list.

    If base_union (shapely geometry of the pre-noise backbone) is given, each
    street is tagged is_noise=True when its constituent geometry is NOT on the
    backbone (i.e. it was added by add_noise_streets).
    """
    from decompose_city import build_dual_from_gdf, decompose
    import networkx as nx

    gdf = segments_to_gdf(segs)
    G_dual = build_dual_from_gdf(gdf, T_BUFFER, already_projected=True,
                                 simplify_roundabout=False)
    res = decompose(G_dual, A_THRESHOLD)

    if base_union is not None:
        # rebuild component membership exactly as decompose() does
        filt = G_dual.copy()
        filt.remove_edges_from([(u, v) for u, v, a in G_dual.edges(data=True)
                                if a['new_angle'] > A_THRESHOLD])
        comps = list(nx.connected_components(filt))
        eps = MIN_LENGTH * 0.01
        is_noise = np.zeros(len(comps), dtype=bool)
        for i, comp in enumerate(comps):
            # noise streets never merge with backbone strokes (angle-protected),
            # so a single off-backbone segment marks the whole street
            for n in comp:
                g = G_dual.nodes[n]['geometry']
                if g.centroid.distance(base_union) > eps:
                    is_noise[i] = True
                    break
        res['is_noise'] = is_noise
    return res


def main():
    from space_filling import run_model, add_noise_streets
    from shapely.geometry import LineString
    from shapely.ops import unary_union

    os.makedirs(OUT_DIR, exist_ok=True)
    for model in (2, 1):
        for seed in (42,):
            t0 = time.time()
            base = run_model(model=model, min_length=MIN_LENGTH,
                             max_iter=500_000, seed=seed, verbose=True)
            print(f'model {model} seed {seed}: {len(base)} segments '
                  f'({time.time()-t0:.0f}s)', flush=True)
            np.save(os.path.join(OUT_DIR, f'segs_m{model}_s{seed}.npy'),
                    np.array(base))
            base_union = unary_union([LineString([a, b]) for a, b in base])

            for nf in (0.0, 0.5, 1.0):
                t0 = time.time()
                if nf == 0.0:
                    segs = base
                else:
                    segs = add_noise_streets(base, noise_frac=nf,
                                             seed=seed ^ 0xBADC0FFE,
                                             min_length=MIN_LENGTH)
                res = decompose_segments(segs,
                                         base_union if nf > 0 else None)
                tag = f'm{model}_s{seed}_nf{nf:g}'
                np.savez(os.path.join(OUT_DIR, f'decomp_{tag}.npz'), **res)
                kmax = res['k'].max() if len(res['k']) else 0
                print(f'  {tag}: {len(res["k"])} streets, kmax={kmax} '
                      f'({time.time()-t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
