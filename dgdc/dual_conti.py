import numpy as np
import math
import networkx as nx
from shapely.geometry import Point, MultiLineString
import geopandas as gpd
import pandas as pd
from shapely.ops import linemerge
import momepy
import osmnx as ox
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def direction(line):
    """Return the (dx, dy) direction vector of a LineString."""
    x, y = line.xy
    return x[-1] - x[0], y[-1] - y[0]

def delta_angle(line1, line2):
    v1 = direction(line1)
    v2 = direction(line2)

    dot = np.dot(v1, v2)
    norm = np.linalg.norm(v1) * np.linalg.norm(v2)

    if norm == 0: # avoids division by zero
        return 0.0

    # acute angle 0-90
    cos_theta = np.clip(dot / norm, -1.0, 1.0)
    return math.degrees(math.acos(abs(cos_theta)))

def new_angles(G, touch_buffer):
    """Compute local angles at touch points and set 'new_angle' on each dual edge."""
    for u, v in G.edges():
        a = G.nodes[u]['geometry']
        b = G.nodes[v]['geometry']
        touch_point = a.intersection(b)

        if isinstance(touch_point, Point):
            touch_point_b = touch_point.buffer(touch_buffer)
            cut_a = a.intersection(touch_point_b)
            cut_b = b.intersection(touch_point_b)
            cut_b = check_string(cut_b, touch_point)
            cut_a = check_string(cut_a, touch_point)
            angle = delta_angle(cut_a, cut_b)
        else: # not a point (e.g parallel lines)
            angle = delta_angle(a, b)

        G[u][v]['new_angle'] = angle

    return G

def check_string(l, p):
    """Return the segment of a MultiLineString that touches point p, or l if already a LineString."""
    if isinstance(l, MultiLineString):
        for i in l.geoms:
            if i.touches(p):
                return i
    else: 
        return l
    
def combine(elements):
    """Flatten a mixed list of ints and lists into a single list."""
    result_list = []
    for item in elements:
        if isinstance(item, int):
            result_list.append(item)
        elif isinstance(item, list):
            result_list.extend(item)
    return result_list


def clean_chains(G_primal):
    """Iteratively merge degree-2 nodes by collapsing them into a single edge."""
    while True:
        nodes_to_remove = []
        
        for node in list(G_primal.nodes()):
            neighbors = list(G_primal.neighbors(node))
            # len(neihbors) is somehow not redundant?
            if  G_primal.degree(node) == 2 and len(neighbors) == 2:
            
                # merge lines and ids
                edge_data_list = list(G_primal.edges(node, data=True))
                my_id = combine([i[2]['id'] for i in edge_data_list])
                lines = [i[2]['geometry'] for i in edge_data_list]
                multi_line = linemerge(MultiLineString(lines))

                # add new edge and remove node
                G_primal.add_edge(neighbors[1], neighbors[0], 
                                id=my_id,
                                geometry=multi_line,
                                new_edge=True)
                nodes_to_remove.append(node)
                
        if not nodes_to_remove: # if empty
            break
        G_primal.remove_nodes_from(nodes_to_remove)

    return G_primal


def split_until_degree_2(G, attr):
    """Remove edges from high-degree nodes (keeping the 2 most similar) until all nodes have degree ≤ 2."""
    G = G.copy()

    changed = True
    while changed:
        changed = False

        for n in list(G.nodes()):
            edges = list(G.edges(n, data=True))
            if len(edges) <= 2:
                continue

            # keep 2 smallest attr (most similar)
            edges_sorted = sorted(edges, key=lambda e: e[2][attr])
            to_remove = edges_sorted[2:]

            for u, v, _ in to_remove:
                if G.has_edge(u, v):
                    G.remove_edge(u, v)
                    changed = True
    return [G.subgraph(c).copy() for c in nx.connected_components(G)]


def merged_G_angle(H, thresh, attr, enforce_degree2):
    """Merge dual-graph nodes with similar angles (within thresh) into single nodes."""
    filtered_H = H.copy()

    # Create components by removing edges with non similar angle
    filtered_H.remove_edges_from([(u, v) for u, v, a in H.edges(data=True) if a[attr] > thresh])
    components = nx.connected_components(filtered_H)
    geometries = nx.get_node_attributes(H, "geometry")
    mapping = {}
    geom_map = {}
    uid_map = {}
    edge_uids = nx.get_node_attributes(H, "edgeUID")

    for comp_nodes in components:
        comp = filtered_H.subgraph(comp_nodes).copy()

        # Optional splitting step
        if enforce_degree2:
            sub_comps = split_until_degree_2(comp, attr)
        else:
            sub_comps = [comp]

        for sub in sub_comps:
            nodes = list(sub.nodes())
            if not nodes:
                continue

            mean_node = tuple(np.mean(np.array(nodes), axis=0))
            for n in nodes:
                mapping[n] = mean_node

            lines = [geometries[n] for n in nodes if n in geometries]
            if lines:
                geom_map[mean_node] = linemerge(MultiLineString(lines))

            uids = []
            for n in nodes:
                val = edge_uids.get(n)
                if val is None:
                    pass
                elif isinstance(val, list):
                    uids.extend(val)
                else:
                    uids.append(val)
            uid_map[mean_node] = uids

    merged_H = nx.relabel_nodes(H, mapping)
    merged_H.remove_edges_from(nx.selfloop_edges(merged_H))  # add this
    nx.set_node_attributes(merged_H, geom_map, "geometry")
    nx.set_node_attributes(merged_H, uid_map, "edgeUID")

    return merged_H

def get_dual_dir_con(t_buffer, a_threshold, data, simplify_roundabout, enforce_degree2):
    """Build the dual continuity graph from an osmnx graph.

    Args:
        t_buffer: buffer radius (meters) around touch points for local angle computation.
        a_threshold: max angle (degrees) for two streets to be considered continuous.
        data: osmnx graph.
        simplify_roundabout: if True, simplify roundabout with momepy/coins
        enforce_degree2: if True, split merged components so no node exceeds degree 2.

    Returns:
        gdf_merged: GeoDataFrame of merged dual nodes with degree and length.
        H: merged dual graph (networkx).
        shape_exploded_df: exploded primal edge GeoDataFrame.
        lines: primal edge GeoDataFrame after chain cleaning.
    """

    if not hasattr(data, "nodes"):
        raise TypeError("data must be an osmnx graph")

    shape_df = ox.graph_to_gdfs(ox.convert.to_undirected(data), nodes=False)
    # shape_df.crs = "epsg:4326" old version, new version below
    shape_df = shape_df.set_crs("epsg:4326")
    estimated_crs = shape_df.estimate_utm_crs()
    shape_df = shape_df.to_crs(estimated_crs)
    
    if simplify_roundabout:
        try:
            shape_df = momepy.roundabout_simplification(shape_df)
            # print('roundabout simplification applied')
        except Exception as e:
            None # assuming there are no roundabouts
            # print(f'roundabout simplification failed: {e}')
        
    if 'edgeUID' not in shape_df.columns:
        shape_df['edgeUID'] = shape_df.index

    # explodes the geometry
    shape_df = shape_df.reset_index().explode('geometry')
    u = shape_df.union_all()
    # i = u.intersection(u) # old version
    i = u
    out = gpd.GeoDataFrame(geometry=gpd.GeoSeries(i, crs=shape_df.crs).explode()).reset_index(drop=True)
    shape_exploded_df = out.sjoin(shape_df[['osmid', 'geometry','edgeUID']], how="left", predicate="intersects")
    shape_exploded_df = shape_exploded_df.drop_duplicates(subset=['geometry'])
    shape_exploded_df = shape_exploded_df.reset_index(drop=True)
    shape_exploded_df['id'] = shape_exploded_df.index

    # converts to primal graph and merges chains 
    G_primal = momepy.gdf_to_nx(shape_exploded_df, approach="primal")
    G_primal = clean_chains(G_primal)
    _, lines = momepy.nx_to_gdf(G_primal)

    # compute angles 
    G_dual = momepy.gdf_to_nx(lines , approach='dual', multigraph=False, angles=False)
    G_dual = new_angles(G_dual,touch_buffer=t_buffer) 

    # merges
    H = merged_G_angle(G_dual,thresh=a_threshold,attr='new_angle',enforce_degree2=enforce_degree2)

    # create dataframe
    df_nodes = pd.DataFrame.from_dict(dict(H.nodes(data=True)), orient='index')
    gdf_merged = gpd.GeoDataFrame(df_nodes, geometry='geometry')
    
    # gdf_merged['degree']=np.array([d for n, d in H.degree()]) # old version, new version below
    degree_map = dict(H.degree())
    gdf_merged['degree'] = gdf_merged.index.map(degree_map)
    gdf_merged['degree_log'] = gdf_merged.degree.apply(lambda x: np.log10(x) if x > 0 else 0)
    gdf_merged['length'] = gdf_merged.geometry.length
    gdf_merged.crs = estimated_crs

    return gdf_merged, H, shape_exploded_df, lines