# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**fractality** (formally "DualGraphDirCon") is a street network morphology research project studying urban street patterns through dual graph continuity analysis and generative space-filling models.

## Setup & Commands

```bash
# Install package in development mode
pip install -e .

# Run city degree sequence analysis (96 cities, 20 parallel processes)
python cities/regime.py

# Run space-filling model (generates model1.shp and model2.shp)
cd space_filling_model && python space_filling.py

# Real-area space-filling model (circular domain sized to a cached city's area)
cd space_filling_model && python real_space_filling.py

# Pre-download OSM graphs to data/city_graphs/ (cached as .graphml)
python cities/download_cities.py

# Grid search over (t_buffer, a_threshold); writes .npz of degree/length arrays
python cities/regime_grid_search.py

# Sierpinski carpet street network + DGDC analysis (writes GeoPackages)
python sierpinski_carpet_model/sierpinski_carpet.py

# Per-intersection street angles for each cached city (writes data/city_angles/<city>.out)
python city_angles/compute_angles.py --t_buffer 50 --processes 20

# COINS-stroke dual-graph degree distribution per city (baseline vs DGDC)
python COINS_comparison/coins_degrees.py --a_threshold 20 --processes 20

# Per-edge compass bearings per city, two methods (writes data/city_orientations/<city>.npz)
python city_orientation/compute_orientations.py --processes 20

# Match Sensible-DTU mobility points to nearest Copenhagen street segment
python mobility/match_sensible.py
```

No automated test suite or linter is configured. Research is done via Jupyter notebooks in `notebooks_lasse/`.

## Architecture

### `dgdc/` — Dual Graph Directional Continuity Package

Core algorithm in `dgdc/dual_conti.py`. Main entry point:

```python
from dgdc import get_dual_dir_con
gdf_merged, H, shape_exploded_df, lines = get_dual_dir_con(
    t_buffer=10,             # buffer radius around intersection touch points for angle computation
    a_threshold=20,          # max angle (degrees) to consider streets as continuous
    data=G,                  # osmnx graph
    simplify_roundabout=True, # collapse roundabouts to nodes via momepy before analysis
    enforce_degree2=False
)
```

**Pipeline:**
1. Fetch OSM street network via `osmnx`
2. `clean_chains()` — iteratively merges degree-2 nodes into single edges
3. `new_angles()` — computes local directional angles at intersections using touch-point buffers
4. `merged_G_angle()` — merges dual-graph nodes with similar angles (below `a_threshold`)
5. Optionally `split_until_degree_2()` — enforces max degree 2 by removing dissimilar edges

### `space_filling_model/space_filling.py` — Generative Street Model

Iterative subdivision model that produces fractal-like street networks:

- **Model 1 (basic):** Select segment proportional to length, split at midpoint, connect to nearest visible intersection
- **Model 2 (biased):** Same as Model 1 with constraints: polygon area `A(r) > 0.05 * exp(-1/r)` and all intersection angles > π/4

Uses `cKDTree` for nearest-neighbor lookups and NumPy-vectorized geometry checks for performance. Dynamic array allocation with doubling capacity.

```python
from space_filling_model.space_filling import run_model, get_geodataframe
segs = run_model(model=2, seed=42, min_angle=np.pi/4)
gdf = get_geodataframe(model=1, seed=42)
```

### `space_filling_model/real_space_filling.py` — Real-Area Space-Filling Model

Wraps `space_filling.run_model` but runs it at real city size on a **circular**
domain: the unit-square is replaced by a circle whose area matches the
convex-hull area of a cached OSM city graph (pass `city=`) or an explicit
`area_m2=`. Init is `"spokes"` (boundary + N lines from the centre). Geometry is
in metres throughout (defaults `min_length=50` m, `min_angle=π/4`). The returned
GeoDataFrame is always georeferenced: the city's local UTM CRS when a `city` is
given, an explicit `crs=` if passed, otherwise `EPSG:32632` (UTM 32N). No water
masking — the domain is a plain circle.

```python
from space_filling_model.real_space_filling import run_circular_model, to_geodataframe
res = run_circular_model(area_m2=np.pi * 500**2, model=2, min_length=50.0, seed=42)
gdf = to_geodataframe(res, model_id=2)
```

### `sierpinski_carpet_model/sierpinski_carpet.py` — Synthetic Fractal Network

Builds a Sierpinski carpet street network (recursive 3×3 subdivision with centre-cell removal; 8^k surviving cells at level `k`) and runs DGDC on it as a ground-truth fractal benchmark.

```python
from sierpinski_carpet_model.sierpinski_carpet import run_dgdc_on_carpet, get_carpet_geodataframe
gdf = get_carpet_geodataframe(level=3)
gdf_carpet, gdf_merged, H = run_dgdc_on_carpet(level=3, a_threshold=20)
```

### `cities/` — City-Scale Analysis

- `cities.txt` — list of ~96 target cities.
- `get_city_list.py` — helper for building/refreshing the target-city list.
- `download_cities.py` — pre-fetches OSM graphs via `osmnx` and caches them as `data/city_graphs/<city>.graphml` (use `load_city(city)` to read back).
- `regime.py` — runs DGDC on all cached cities in parallel (20 processes); `(t_buffer, a_threshold)` are CLI args (default `10`/`20`); writes `data/city_degrees/t<t>_a<a>/<city>.out`.
- `regime_grid_search.py` — sweeps `(t_buffer, a_threshold)` combinations and writes `.npz` files containing degree and length arrays per city.

### `city_angles/compute_angles.py` — Per-Intersection Angle Dump

Re-runs the DGDC pipeline up to `new_angles` for every cached `data/city_graphs/<city>.graphml`, then writes the flat list of dual-graph edge angles (acute, 0–90°) to `data/city_angles/<city>.out` via `np.savetxt`. Runs in a `multiprocessing.Pool` (default 20 procs); existing outputs are skipped.

### `COINS_comparison/coins_degrees.py` — COINS Baseline

Builds a dual graph from `momepy.COINS` strokes (each stroke = a node; strokes
sharing a primal endpoint are connected) and writes the `(degree, length)`
arrays per city to `data/city_degrees_coins/a<a>/<city>.npz`. Provides a
continuity-stroke baseline to compare against DGDC degree distributions. COINS's
`angle_threshold` is an interior angle, so DGDC's max-deflection `a_threshold` is
passed as `180 - a_threshold`. `compare_degrees.ipynb` does the comparison.

### `city_orientation/compute_orientations.py` — Street Bearings

Computes per-edge compass bearings for each cached city graph (Boeing 2019),
offline. Method 1 uses simplified-topology endpoint-to-endpoint bearings (plus
reciprocal); Method 2 uses the raw straight-line pieces of each edge's
`LineString` geometry, each weighted by its metric length. Writes
`data/city_orientations/<city>.npz` (`method1_bearings`, `method2_bearings`,
`method2_weights`). `orientations.ipynb` analyses the output.

### `mobility/match_sensible.py` — Trajectory-to-Network Matching

Matches Sensible-DTU mobility points (Greater Copenhagen) to the nearest DGDC
street segment via a `shapely.STRtree` nearest-neighbour query, attaching each
point's segment id, dual-graph `degree`, and `osmid`. Splits points into
moving (`label == -1`) and stop (infostop) sets and writes
`moving_segments.csv` / `stops_segments.csv`. `sensible.ipynb` explores the
result.

### Output directories

`data/city_graphs/` holds cached `.graphml` OSM graphs; `data/city_degrees/` holds DGDC degree-sequence outputs from `regime.py`; `data/city_degrees_coins/` holds the COINS baseline degree/length `.npz` files; `data/city_angles/` holds per-edge angle arrays from `city_angles/compute_angles.py`; `data/city_orientations/` holds per-city bearing `.npz` files; `data/city_external/` holds `city_info.csv` (contextual per-city attributes — population, founding era, street pattern, climate, etc., documented in its own README); `data/traj/` holds additional trajectory data. `cache/` holds JSON osmnx HTTP response cache.

## Key Dependencies

| Package | Purpose |
|---------|---------|
| `osmnx` | Fetch OSM street networks |
| `momepy` | Urban morphology (roundabout simplification, primal/dual graph conversion) |
| `networkx` | Graph construction and analysis |
| `shapely` | Geometric operations (LineString intersections, angles) |
| `scipy.spatial.cKDTree` | Spatial nearest-neighbor lookups |
| `geopandas` | Spatial DataFrames for results |
