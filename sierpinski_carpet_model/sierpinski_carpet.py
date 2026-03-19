#!/usr/bin/env python3
"""
Sierpinski carpet street network model.

Generates a Sierpinski carpet as a planar street network at a given recursion
level and runs the dgdc dual-graph directional continuity algorithm on it.

The carpet is built by the standard recursive removal: start with the unit
square, divide into a 3×3 grid, remove the centre cell, and repeat for each of
the 8 remaining cells.  The boundaries of the surviving cells form the street
network (horizontal and vertical segments).

Usage
-----
    from sierpinski_carpet import run_dgdc_on_carpet, get_carpet_geodataframe

    gdf_carpet, gdf_merged, H = run_dgdc_on_carpet(level=3, a_threshold=20)
"""

import os
import sys
import numpy as np
from shapely.geometry import LineString
from collections import defaultdict
import geopandas as gpd

# ---------------------------------------------------------------------------
# Sierpinski carpet geometry
# ---------------------------------------------------------------------------

def get_carpet_cells(level):
    """
    Return the set of (i, j) integer cell indices present in the Sierpinski
    carpet at the given recursion level.  The grid has 3^level × 3^level cells.
    """
    if level == 0:
        return {(0, 0)}
    prev = get_carpet_cells(level - 1)
    result = set()
    for (pi, pj) in prev:
        for di in range(3):
            for dj in range(3):
                if di == 1 and dj == 1:
                    continue  # remove centre
                result.add((pi * 3 + di, pj * 3 + dj))
    return result


def get_carpet_segments(level):
    """
    Generate the unique line segments that form the boundary network of the
    Sierpinski carpet at the given recursion level.

    Returns a list of ((x0, y0), (x1, y1)) tuples in the unit square [0,1]².
    Each segment is a side of a surviving cell; shared sides between adjacent
    cells appear only once.
    """
    cells = get_carpet_cells(level)
    n = 3 ** level
    s = 1.0 / n

    segments = set()
    for (i, j) in cells:
        x0 = round(i * s, 12)
        y0 = round(j * s, 12)
        x1 = round((i + 1) * s, 12)
        y1 = round((j + 1) * s, 12)
        for edge in [
            ((x0, y0), (x1, y0)),  # bottom
            ((x0, y1), (x1, y1)),  # top
            ((x0, y0), (x0, y1)),  # left
            ((x1, y0), (x1, y1)),  # right
        ]:
            segments.add(tuple(sorted(edge)))

    return list(segments)


def to_geodataframe(segments):
    """Convert a list of ((x0,y0),(x1,y1)) tuples to a GeoDataFrame."""
    records = [
        {'geometry': LineString([a, b]), 'osmid': i}
        for i, (a, b) in enumerate(segments)
    ]
    return gpd.GeoDataFrame(records, crs='EPSG:4326')


def get_carpet_geodataframe(level=3):
    """
    Build and return a GeoDataFrame of Sierpinski carpet street segments at
    the given recursion level.

    Parameters
    ----------
    level : int
        Recursion depth.  Level k produces 8^k surviving cells and at most
        4 × 8^k unique segment edges (many are shared).
        Practical range: 1–4.
    """
    segs = get_carpet_segments(level)
    return to_geodataframe(segs)


# ---------------------------------------------------------------------------
# Degree helpers
# ---------------------------------------------------------------------------

def compute_degrees(segments):
    """Return a dict mapping each rounded node key to its degree."""
    deg = defaultdict(int)
    for a, b in segments:
        ka = (round(a[0], 10), round(a[1], 10))
        kb = (round(b[0], 10), round(b[1], 10))
        deg[ka] += 1
        deg[kb] += 1
    return dict(deg)


# ---------------------------------------------------------------------------
# Run dgdc on the carpet
# ---------------------------------------------------------------------------

def run_dgdc_on_carpet(level=3, t_buffer=1, a_threshold=20,
                       enforce_degree2=False):
    """
    Generate the Sierpinski carpet network and run the dual-graph directional
    continuity algorithm.

    Parameters
    ----------
    level           : recursion depth (2 = 64 segments, 3 = 512 segments)
    t_buffer        : touch-point buffer radius (in projected metres) for the
                      local angle computation inside dgdc
    a_threshold     : max angle (degrees) for two streets to count as
                      directionally continuous
    enforce_degree2 : if True, split merged components so no dual node exceeds
                      degree 2

    Returns
    -------
    gdf_carpet  : GeoDataFrame of primal carpet segments
    gdf_merged  : GeoDataFrame of merged dual nodes (degree, length, …)
    H           : merged dual graph (networkx.Graph)
    """
    import momepy

    # Ensure the repo root is importable regardless of working directory
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)

    from dgdc import get_dual_dir_con

    print(f"Building carpet at level {level}…")
    gdf_carpet = get_carpet_geodataframe(level)
    print(f"  {len(gdf_carpet)} primal segments")

    G = momepy.gdf_to_nx(gdf_carpet, approach='primal')

    print("Running dgdc…")
    gdf_merged, H, _, _ = get_dual_dir_con(
        t_buffer=t_buffer,
        a_threshold=a_threshold,
        data=G,
        enforce_degree2=enforce_degree2,
        simplify_roundabout=False,
    )
    print(f"  {H.number_of_nodes()} dual nodes, {H.number_of_edges()} dual edges")
    return gdf_carpet, gdf_merged, H


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("=== Sierpinski Carpet Street Network Model ===\n")

    for level in [2, 3]:
        gdf_carpet, gdf_merged, H = run_dgdc_on_carpet(level=level)
        out = f'carpet_level{level}_dual.gpkg'
        gdf_merged.to_file(out, driver='GPKG')
        print(f"  Saved dual nodes → {out}\n")
