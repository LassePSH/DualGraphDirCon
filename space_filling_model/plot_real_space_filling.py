import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm
import geopandas as gpd
import momepy
from tqdm import tqdm
from shapely.ops import unary_union
from real_space_filling import run_circular_model, to_geodataframe

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from dgdc import get_dual_dir_con

print("Running real space-filling model and plotting results...")

# p_stub injects random dead-end stubs to populate the low-degree (k<=4)
# shoulder of the dual-graph degree distribution. Sweep it (e.g. 0.0, 0.3,
# 0.6, 0.9) and re-run to see the shoulder rise; lower snap_tol (e.g. 1.0) for
# more degree-1 cul-de-sacs instead of degree-2 T-ins.
P_STUB = 0.3
res = run_circular_model(area_m2=np.pi*5000**2, area_coeff = 1_000, model=2, min_length=50.0, seed=42, p_stub=P_STUB)
gdf = to_geodataframe(res, model_id=2)

# Per-segment degree: the higher of its two endpoint degrees. This highlights


print("Running DGDC on the model output to obtain the dual-graph degree sequence...")
gdf_dgdc = gdf[['geometry']].copy()
gdf_dgdc['osmid'] = np.arange(len(gdf_dgdc))
gdf_dgdc = gdf_dgdc.to_crs(4326)

G = momepy.gdf_to_nx(gdf_dgdc, approach='primal')
gdf_merged, H, _, _ = get_dual_dir_con(t_buffer=1, a_threshold=20, data=G,
                                       enforce_degree2=False,
                                       simplify_roundabout=False)


OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')
os.makedirs(OUT_DIR, exist_ok=True)


# --------------------------------------------------------------------------- #
# 1. Raw lines                                                                #
# --------------------------------------------------------------------------- #
# plot the raw lines. make sure the lines are visible by setting a fitting
# linewidth and alpha.
fig, ax = plt.subplots(figsize=(10, 10))
gdf.plot(ax=ax, color='black', linewidth=0.4, alpha=0.7)
ax.set_aspect('equal')
ax.set_axis_off()
ax.set_title(f'Real space-filling model (model {gdf["model"].iloc[0]}) — 'f'{len(gdf)} segments')
fig.tight_layout()
raw_path = os.path.join(OUT_DIR, 'real_space_filling_raw_2.png')
fig.savefig(raw_path, dpi=200, bbox_inches='tight')
print(f'wrote {raw_path}')


# --------------------------------------------------------------------------- #
# 2. Segments coloured by degree                                              #
# --------------------------------------------------------------------------- #
# plot the segments colored by degree. Bin the degree into 10 bins and use a
# colormap to color the segments accordingly.
# plot the low degree first and the high degree last to make sure the high
# degree segments are visible.
N_BINS = 10
cmap = plt.get_cmap('viridis', N_BINS)
dmin, dmax = gdf_merged.degree_log.min(), gdf_merged.degree_log.max()
# Integer-aware bin edges; fall back to a linear split if the degree range is
# narrower than the requested number of bins.
edges = np.linspace(dmin, dmax + 1, N_BINS + 1)
norm = BoundaryNorm(edges, cmap.N)
# Draw low degree first, high degree last so busy segments sit on top.
gdf_merged = gdf_merged.sort_values('degree_log', kind='stable')

fig, ax = plt.subplots(figsize=(10, 10))
gdf_merged.plot(
    ax=ax,
    column='degree_log',
    cmap=cmap,
    norm=norm,
    linewidth=0.6,
    alpha=0.9,
    legend=False,
)
ax.set_aspect('equal')
ax.set_axis_off()
ax.set_title('Segments coloured by endpoint degree')

sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('degree')

fig.tight_layout()
deg_path = os.path.join(OUT_DIR, 'real_space_filling_degree_2.png')
fig.savefig(deg_path, dpi=100, bbox_inches='tight')
print(f'wrote {deg_path}')


# --------------------------------------------------------------------------- #
# 3. Dual-graph degree distribution for the segments above                    #
# --------------------------------------------------------------------------- #
# Run DGDC on the model output produced for the plots above to obtain the
# dual-graph (directional-continuity) degree sequence, then plot p(k) vs k.
#
# `get_dual_dir_con` forces the input graph to EPSG:4326 and then estimates a
# UTM CRS, so it expects geographic (lat/lon) coordinates. The model `gdf` is
# in UTM metres, so we reproject it back to 4326 before handing it over;
# t_buffer is then a meaningful metre buffer in the re-estimated UTM CRS.
# `simplify_roundabout=False` matches the synthetic-network usage in the
# Sierpinski carpet model (there are no OSM roundabouts to collapse here).
print(f'DGDC: {H.number_of_nodes()} dual nodes, {H.number_of_edges()} dual edges')

d = np.asarray(gdf_merged['degree'], dtype=float)
d = d[d > 0]

fig, ax = plt.subplots(figsize=(7, 6))

# Log-spaced bins, one per distinct degree value (at least 2 edges).
n_bins = max(len(np.unique(d)), 2)
y, x = np.histogram(d, bins=np.logspace(0, np.log10(d.max()), n_bins), density=True)
x = x[:-1]
x_nn = x[y != 0]
y_nn = y[y != 0]
ax.plot(x_nn, y_nn, marker='o', ms=4, alpha=1, zorder=2, color='k')

ax.yaxis.set_minor_formatter(plt.NullFormatter())

# Hide top/right spines (Nature style)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only bottom/left axes
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# Major ticks (not too many)
ax.xaxis.set_major_locator(plt.LogLocator(base=10.0))
ax.yaxis.set_major_locator(plt.LogLocator(base=10.0))

# Minor ticks
ax.xaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=np.arange(1, 10)))
ax.yaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=np.arange(1, 10)))

# Appearance
ax.tick_params(which='major', length=6, width=1.2)
ax.tick_params(which='minor', length=3, width=1.0)

ax.axvline(4, color='k', linestyle='dotted', label=r'$\phi = $' + str(4), zorder=-1)
plt.legend()

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_title("Degree Distribution", loc='left', pad=12)
ax.set_xlabel("k")
ax.set_ylabel("p(k)")

fig.tight_layout()
dist_path = os.path.join(OUT_DIR, 'real_space_filling_degree_dist_2.png')
fig.savefig(dist_path, dpi=100, bbox_inches='tight')
print(f'wrote {dist_path}')

plt.show()
