"""Run DGDC on the full Denmark OSM extract and plot the dual-graph degrees.

The input is a raw OpenStreetMap ``.osm.pbf`` country extract. ``pyrosm`` reads
it and converts the street network into an osmnx-compatible graph, which is then
fed to ``dgdc.get_dual_dir_con``.

Because both steps are expensive on a whole-country dataset, intermediate
results are cached:

* the converted osmnx graph is saved as ``.graphml`` (so the slow PBF scan only
  happens once), and
* the DGDC output GeoDataFrame is saved as GeoParquet (so the plot can be
  re-tuned without re-running DGDC).

Examples
--------
    # full country (uses/creates caches next to this script)
    python dgdc_denmark.py

    # quick test on a small bounding box (minx miny maxx maxy, lon/lat)
    python dgdc_denmark.py --bbox 12.56 55.67 12.60 55.69 --tag cph

    # only redo the figure from a cached DGDC result
    python dgdc_denmark.py --reuse
"""

import os
import sys
import time
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless: we only write files
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm

import geopandas as gpd
import osmnx as ox
from pyrosm import OSM

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from dgdc import get_dual_dir_con

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PBF = "/work/lpsha/data/OSM_denmark/denmark-latest.osm.pbf"


def pbf_to_osmnx_graph(pbf_path, network_type="all", bounding_box=None,
                       cache_path=None, reuse=True):
    """Read a .osm.pbf and return an osmnx-compatible MultiDiGraph (EPSG:4326).

    The converted graph is cached as GraphML at ``cache_path`` so the (slow) PBF
    parse is only done once. Pass ``reuse=False`` to force a rebuild.
    """
    if cache_path and reuse and os.path.exists(cache_path):
        print(f"loading cached graph from {cache_path}")
        return ox.load_graphml(cache_path)

    if not os.path.exists(pbf_path):
        raise FileNotFoundError(f"PBF not found: {pbf_path}")

    print(f"reading network ({network_type}) from {pbf_path} ...")
    t = time.time()
    osm = OSM(pbf_path, bounding_box=bounding_box)
    nodes, edges = osm.get_network(network_type=network_type, nodes=True)
    print(f"  parsed {len(edges)} edges in {time.time() - t:.0f}s; building graph ...")

    G = osm.to_graph(nodes, edges, graph_type="networkx", osmnx_compatible=True)
    print(f"  graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
          f"crs={G.graph.get('crs')}")

    if cache_path:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        print(f"  caching graph -> {cache_path}")
        ox.save_graphml(G, cache_path)
    return G


def plot_degrees(gdf_merged, out_path, n_bins=10, figsize=20.0, dpi=300,
                 linewidth=0.3):
    """Plot dual-graph segments coloured by endpoint degree (log scale).

    Tuned for very large datasets: the line geometry is rasterised (so the PNG
    does not embed millions of vector paths) and antialiasing is disabled.
    """
    cmap = plt.get_cmap("viridis", n_bins)
    dmin, dmax = gdf_merged.degree_log.min(), gdf_merged.degree_log.max()
    edges = np.linspace(dmin, dmax + 1, n_bins + 1)
    norm = BoundaryNorm(edges, cmap.N)

    # draw low-degree first so high-degree (rare, important) ends up on top
    gdf_merged = gdf_merged.sort_values("degree_log", kind="stable")

    fig, ax = plt.subplots(figsize=(figsize, figsize))
    print(f"rendering {len(gdf_merged)} segments ...")
    t = time.time()
    gdf_merged.plot(
        ax=ax,
        column="degree_log",
        cmap=cmap,
        norm=norm,
        linewidth=linewidth,
        legend=False,
        antialiased=False,  # faster for millions of lines
        rasterized=True,    # keep the saved figure light
    )
    print(f"  plotted in {time.time() - t:.0f}s")

    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("Denmark — segments coloured by dual-graph endpoint degree")

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("degree (log10)")

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pbf", default=DEFAULT_PBF, help="path to the .osm.pbf extract")
    p.add_argument("--network-type", default="all",
                   help="pyrosm network filter (all, driving, walking, ...)")
    p.add_argument("--bbox", type=float, nargs=4, metavar=("MINX", "MINY", "MAXX", "MAXY"),
                   default=None, help="optional lon/lat bounding box to subset (for testing)")
    p.add_argument("--t-buffer", type=float, default=10,
                   help="DGDC touch-point buffer radius in metres")
    p.add_argument("--a-threshold", type=float, default=20.0,
                   help="DGDC max continuity angle in degrees")
    p.add_argument("--tag", default="denmark",
                   help="name stem for cache + output files")
    p.add_argument("--reuse", action="store_true",
                   help="reuse cached graph/DGDC result if present (else rebuild)")
    p.add_argument("--n-bins", type=int, default=10)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--figsize", type=float, default=20.0)
    p.add_argument("--linewidth", type=float, default=0.3)
    args = p.parse_args()

    cache_dir = os.path.join(HERE, "cache")
    graph_cache = os.path.join(cache_dir, f"{args.tag}_{args.network_type}.graphml")
    dgdc_cache = os.path.join(cache_dir, f"{args.tag}_dgdc.parquet")
    out_png = os.path.join(HERE, f"{args.tag}_dgdc.png")

    # 1. DGDC result: reuse cached GeoParquet if asked, else (re)compute.
    if args.reuse and os.path.exists(dgdc_cache):
        print(f"loading cached DGDC result from {dgdc_cache}")
        gdf_merged = gpd.read_parquet(dgdc_cache)
    else:
        G = pbf_to_osmnx_graph(
            args.pbf,
            network_type=args.network_type,
            bounding_box=args.bbox,
            cache_path=graph_cache,
            reuse=args.reuse,
        )
        print("running DGDC ...")
        t = time.time()
        gdf_merged, H, _, _ = get_dual_dir_con(
            t_buffer=args.t_buffer,
            a_threshold=args.a_threshold,
            data=G,
            enforce_degree2=False,
            simplify_roundabout=True,
        )
        print(f"  DGDC produced {len(gdf_merged)} segments in {time.time() - t:.0f}s")
        os.makedirs(cache_dir, exist_ok=True)
        # Keep only what the plot needs: the merged graph also carries object
        # columns (e.g. lists of merged node ids) that Parquet cannot encode.
        slim = gdf_merged[["geometry", "degree", "degree_log", "length"]]
        slim.to_parquet(dgdc_cache)
        gdf_merged = slim
        print(f"  cached DGDC result -> {dgdc_cache}")

    # 2. Plot.
    plot_degrees(
        gdf_merged,
        out_png,
        n_bins=args.n_bins,
        figsize=args.figsize,
        dpi=args.dpi,
        linewidth=args.linewidth,
    )


if __name__ == "__main__":
    main()
