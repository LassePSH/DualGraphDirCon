"""E2: Instrumented DGDC decomposition.

Re-runs the DGDC pipeline (same steps as dgdc.get_dual_dir_con) but keeps track
of which cleaned primal segments compose each merged dual street, and splits
each street's dual degree k into:

  k_end : distinct neighbor streets attached ONLY at the street's terminal
          nodes (its two free ends),
  k_int : distinct neighbor streets attached at interior junctions, i.e.
          junctions the street passes THROUGH.

Also records: n_seg (segments merged), length, number of terminals, number of
dead-end terminals, per-terminal distinct-neighbor counts, and the primal
junction valence distribution.

Hypothesis under test: k_end <= ~4 (2 ends x <=2 distinct streets/end), so any
street with k > 4 must have interior junctions => the transition phi = 4.
"""
import os
import sys
import numpy as np

sys.path.insert(0, '/home/lpsha/s154446/fractality')

GRAPH_DIR = '/home/lpsha/s154446/fractality/data/city_graphs'
OUT_DIR = '/home/lpsha/s154446/fractality/the_answer/results/decomp_t10_a20'


def build_dual_from_gdf(shape_df, t_buffer, already_projected=False,
                        simplify_roundabout=True):
    """Mirror of get_dual_dir_con up to the angle-annotated dual graph."""
    import geopandas as gpd
    import momepy
    from dgdc.dual_conti import clean_chains, new_angles

    if not already_projected:
        shape_df = shape_df.set_crs('epsg:4326')
        estimated_crs = shape_df.estimate_utm_crs()
        shape_df = shape_df.to_crs(estimated_crs)
    if simplify_roundabout:
        try:
            shape_df = momepy.roundabout_simplification(shape_df)
        except Exception:
            pass
    if 'edgeUID' not in shape_df.columns:
        shape_df['edgeUID'] = shape_df.index
    shape_df = shape_df.reset_index().explode('geometry')
    u = shape_df.union_all()
    out = gpd.GeoDataFrame(geometry=gpd.GeoSeries(u, crs=shape_df.crs)
                           .explode()).reset_index(drop=True)
    cols = [c for c in ('osmid', 'geometry', 'edgeUID') if c in shape_df.columns
            or c == 'geometry']
    sed = out.sjoin(shape_df[cols], how='left', predicate='intersects')
    sed = sed.drop_duplicates(subset=['geometry']).reset_index(drop=True)
    sed['id'] = sed.index

    G_primal = momepy.gdf_to_nx(sed, approach='primal')
    G_primal = clean_chains(G_primal)
    _, lines = momepy.nx_to_gdf(G_primal)
    G_dual = momepy.gdf_to_nx(lines, approach='dual', multigraph=False,
                              angles=False)
    G_dual = new_angles(G_dual, touch_buffer=t_buffer)
    return G_dual


def endpoints_of(geom, nd=3):
    """Rounded coordinate tuples of the free ends of a (Multi)LineString."""
    b = geom.boundary
    if b.is_empty:
        return []
    if b.geom_type == 'Point':
        return [(round(b.x, nd), round(b.y, nd))]
    return [(round(p.x, nd), round(p.y, nd)) for p in b.geoms]


def decompose(G_dual, a_threshold):
    """Split each merged street's dual degree into endpoint/interior parts."""
    import networkx as nx

    nodes = list(G_dual.nodes())
    geom = {n: G_dual.nodes[n]['geometry'] for n in nodes}
    ends = {n: endpoints_of(geom[n]) for n in nodes}

    # global incidence: junction -> segments ending there (primal valence)
    inc = {}
    for n in nodes:
        for p in ends[n]:
            inc.setdefault(p, []).append(n)

    # streets = components of the angle-filtered dual graph
    filt = G_dual.copy()
    filt.remove_edges_from([(u, v) for u, v, a in G_dual.edges(data=True)
                            if a['new_angle'] > a_threshold])
    comps = list(nx.connected_components(filt))
    comp_id = {}
    for i, comp in enumerate(comps):
        for n in comp:
            comp_id[n] = i
    n_streets = len(comps)

    # per-street endpoint multiplicity m_S(p)
    mult = [dict() for _ in range(n_streets)]
    for n in nodes:
        ci = comp_id[n]
        for p in ends[n]:
            mult[ci][p] = mult[ci].get(p, 0) + 1

    # neighbor classification via inter-component dual edges
    nbr_interior = [set() for _ in range(n_streets)]   # neighbors at interior
    nbr_endpoint = [set() for _ in range(n_streets)]   # neighbors at terminals
    nbr_unknown = [set() for _ in range(n_streets)]    # no shared endpoint found
    for u, v in G_dual.edges():
        cu, cv = comp_id[u], comp_id[v]
        if cu == cv:
            continue
        shared = set(ends[u]) & set(ends[v])
        if not shared:
            nbr_unknown[cu].add(cv)
            nbr_unknown[cv].add(cu)
            continue
        for p in shared:
            if mult[cu].get(p, 0) >= 2:
                nbr_interior[cu].add(cv)
            else:
                nbr_endpoint[cu].add(cv)
            if mult[cv].get(p, 0) >= 2:
                nbr_interior[cv].add(cu)
            else:
                nbr_endpoint[cv].add(cu)

    # per-street arrays
    k = np.zeros(n_streets, dtype=np.int64)
    k_int = np.zeros(n_streets, dtype=np.int64)
    k_end = np.zeros(n_streets, dtype=np.int64)
    k_unk = np.zeros(n_streets, dtype=np.int64)
    n_seg = np.array([len(c) for c in comps], dtype=np.int64)
    length = np.zeros(n_streets)
    n_term = np.zeros(n_streets, dtype=np.int64)
    n_dead = np.zeros(n_streets, dtype=np.int64)
    per_end_counts = []

    for i, comp in enumerate(comps):
        interior = nbr_interior[i]
        endp = nbr_endpoint[i] - interior
        unk = nbr_unknown[i] - interior - endp
        k_int[i] = len(interior)
        k_end[i] = len(endp)
        k_unk[i] = len(unk)
        k[i] = len(interior | endp | unk)
        length[i] = sum(geom[n].length for n in comp)
        terms = [p for p, m in mult[i].items() if m == 1]
        n_term[i] = len(terms)
        for p in terms:
            others = {comp_id[n] for n in inc[p] if comp_id[n] != i}
            per_end_counts.append(len(others))
            if len(inc[p]) == 1:
                n_dead[i] += 1

    valence = np.array([len(v) for v in inc.values()], dtype=np.int64)
    return dict(k=k, k_int=k_int, k_end=k_end, k_unk=k_unk, n_seg=n_seg,
                length=length, n_term=n_term, n_dead=n_dead,
                per_end_counts=np.array(per_end_counts, dtype=np.int64),
                valence=valence)


def run_city(city, t_buffer=10, a_threshold=20, out_dir=OUT_DIR):
    import time
    import osmnx as ox
    out_path = os.path.join(out_dir, city + '.npz')
    if os.path.exists(out_path):
        print(city, 'already done'); return
    t0 = time.time()
    fp = os.path.join(GRAPH_DIR, city.replace(' ', '_') + '.graphml')
    G = ox.load_graphml(fp)
    shape_df = ox.graph_to_gdfs(ox.convert.to_undirected(G), nodes=False)
    G_dual = build_dual_from_gdf(shape_df, t_buffer)
    res = decompose(G_dual, a_threshold)
    np.savez(out_path, **res)
    print(f'{city}: {len(res["k"])} streets, {time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    from multiprocessing import Pool

    cities = sys.argv[1:] if len(sys.argv) > 1 else ['Luxembourg City, Luxembourg']
    os.makedirs(OUT_DIR, exist_ok=True)
    if len(cities) == 1:
        run_city(cities[0])
    else:
        with Pool(min(len(cities), 14)) as p:
            p.map(run_city, cities)
