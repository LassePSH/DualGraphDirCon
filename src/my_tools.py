import geopandas as gpd
import pickle
import pandas as pd
import numpy as np
import ast
import math
import networkx as nx
from shapely.geometry import Point, Polygon, MultiPoint, LineString, MultiLineString
from tqdm import tqdm
from shapely import geometry, ops
import shapely
import momepy

# OUTDATED 
# see util.py

#### avg angel ####
def get_slope(line):
    l = len(line.xy[0])-1
    x1, y1, x2, y2 = line.xy[0][0], line.xy[1][0], line.xy[0][l], line.xy[1][l]
    if x2 == x1: # straight vertical line
        return np.inf
    return (y2 - y1) / (x2 - x1)

    
def delta_angle(line1, line2):
    slope1 = get_slope(line1)
    slope2 = get_slope(line2)

    if slope1 == np.inf and slope2 == np.inf:
        return 0
    
    elif slope1 == np.inf:
        return 90-math.degrees(abs(math.atan(slope2)))

    elif slope2 == np.inf:
        return 90-math.degrees(abs(math.atan(slope1)))
    
    elif slope1 * slope2 == -1:
        return 90

    else:
        angle_radians = math.atan(abs((slope2 - slope1) / (1 + slope1 * slope2)))
        angle_degrees = math.degrees(angle_radians)

        return angle_degrees


def new_angles(G):
    nodes = G.nodes(data=True)
    for edges in list(G.edges()):
        u, v = edges
        a =  nodes[u]['geometry']
        b =  nodes[v]['geometry']
        touch_point = a.intersection(b)
        
        if type(touch_point)!=Point:
            angle=delta_angle(a,b)
        else:  
            # thresh_cut = 0.00075633427357771 / 10 #mean(length) / N
            # thresh_cut = 10
            thresh_cut = 5 #meters?
            touch_point_b=touch_point.buffer(thresh_cut*2)
            cut_a = a.difference(a.difference(touch_point_b)) 
            cut_b = b.difference(b.difference(touch_point_b))
            cut_b = check_string(cut_b,touch_point)
            cut_a = check_string(cut_a,touch_point)
            angle=delta_angle(cut_a, cut_b)
        
        G[u][v]['new_angle'] = angle

    return G

def check_string(l,p):
    if type(l) == MultiLineString:
        for i in l.geoms:
            if i.touches(p):
                return i
    else: 
        return l


def merged_G_angle(G,tresh):
    H = G.copy() 

    # find components that have a similar angle and merge them
    filtered_H = H.copy()
    edges_to_remove = [(u, v) for u, v, a in H.edges(data=True) if (a['new_angle'] > tresh) or (a['new_angle'] > (180 - tresh))]
    filtered_H.remove_edges_from(edges_to_remove)
    
    S = [H.subgraph(c).copy() for c in nx.connected_components(filtered_H)]

    lengths = nx.get_node_attributes(G, "mm_len")

    mapping = {}
    length_map = {}
    for s in S: # components
        node_list = list(s.nodes())
        mean_node = tuple(np.mean(np.array(node_list), axis=0))
        mapping.update({n: mean_node for n in node_list})
        length_map[mean_node] = np.array([lengths.get(key) for key in node_list]).sum()

    # relabel the nodes in the original graph
    H=nx.relabel_nodes(H, mapping)
    return H, mapping, length_map

def get_degree(graph,id):
    try:
        return graph.degree(id)
    except: return np.nan

def combine(elements):
    result_list = []
    for item in elements:
        if isinstance(item, int):
            result_list.append(item)
        elif isinstance(item, list):
            result_list.extend(item)

    return result_list


def clean_chains(G_primal):
    while True:
        nodes_to_remove = []
        for node in G_primal.nodes():
            neighbors = list(G_primal.neighbors(node))
            degree = G_primal.degree(node)
            if len(neighbors) == 2 and degree == 2:

                osmIDS = combine([i[2]['id'] for i in list(G_primal.edges(node, data=True))])
                lines = [i[2]['geometry'] for i in list(G_primal.edges(node, data=True))]

                multi_line = ops.linemerge(geometry.MultiLineString(lines))

                G_primal.add_edge(neighbors[1], neighbors[0], 
                                id=osmIDS,
                                geometry=multi_line,
                                new_edge=True)

                nodes_to_remove.append(node)
        if not nodes_to_remove: # if empty
            break
        G_primal.remove_nodes_from(nodes_to_remove)

    return G_primal