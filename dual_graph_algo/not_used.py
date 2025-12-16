# old version

def merged_G_angle(G,thresh,attr):
    H = G.copy() 

    # find components that have a similar angle and merge them
    filtered_H = H.copy()
    edges_to_remove = [(u, v) for u, v, a in H.edges(data=True) if (a[attr] > thresh)]  # or (a[attr] > (180 - tresh))
    filtered_H.remove_edges_from(edges_to_remove)

    # Find connected components (groups of nodes to merge)
    components = [H.subgraph(c).copy() for c in nx.connected_components(filtered_H)]
    geometries = nx.get_node_attributes(G, "geometry")
    
    mapping = {}
    geom_map = {}
    
    for comp in components:
        nodes = list(comp.nodes())
        
        # new node id mean node coordinate
        mean_node = tuple(np.mean(np.array(nodes), axis=0))
        mapping.update({n: mean_node for n in nodes})
    
        # Merge geometries (LineStrings → MultiLineString → merged line)
        lines = [geometries[n] for n in nodes if n in geometries]
        if lines:
            merged_geom = linemerge(MultiLineString(lines))
            geom_map[mean_node] = merged_geom
    
    # Relabel (merge) nodes in the original graph
    merged_H = nx.relabel_nodes(H, mapping)
    nx.set_node_attributes(merged_H, geom_map, "geometry")
    
    return merged_H
