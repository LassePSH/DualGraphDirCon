# Dual Graph Directional Continuity — Pseudocode

```
FUNCTION get_dual_dir_con(streets, angle_threshold):

    streets = reproject(streets) → UTM
    streets = simplify(streets)

    FOR each degree-2 node in streets:
        merge with neighbors into one segment

    dual_graph = convert_to_dual_graph(streets) → nodes are segments, edges are intersections

    FOR each edge in dual_graph:
        angle[edge] = measure_angle(segment_A, segment_B, at intersection)

    FOR each group of segments where angle < angle_threshold:
        collapse group → single "street axis" node

    RETURN dual_graph with degree and length per axis
```
