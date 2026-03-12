# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**fractality** (formally "DualGraphDirCon") is a street network morphology research project studying urban street patterns through dual graph continuity analysis and generative space-filling models.

## Setup & Commands

```bash
# Install package in development mode
pip install -e .

# Run city degree sequence analysis (84 cities, 20 parallel processes)
python cities/regime.py

# Run space-filling model (generates model1.shp and model2.shp)
cd space_filling_model && python space_filling.py
```

No automated test suite or linter is configured. Research is done via Jupyter notebooks in `notebooks_lasse/`.

## Architecture

### `dgdc/` — Dual Graph Directional Continuity Package

Core algorithm in `dgdc/dual_conti.py`. Main entry point:

```python
from dgdc import get_dual_dir_con
gdf_merged, H, shape_exploded_df, lines = get_dual_dir_con(
    t_buffer=10,        # buffer radius around intersection touch points for angle computation
    a_threshold=20,     # max angle (degrees) to consider streets as continuous
    data=G,             # osmnx graph
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

### `cities/regime.py` — City Analysis

Fetches ~84 city street networks via `osmnx`, computes dual graph degree sequences in parallel (20 processes), saves results to `data/city_degrees/t10_a20/<city>.out`.

## Key Dependencies

| Package | Purpose |
|---------|---------|
| `osmnx` | Fetch OSM street networks |
| `momepy` | Urban morphology (roundabout simplification, primal/dual graph conversion) |
| `networkx` | Graph construction and analysis |
| `shapely` | Geometric operations (LineString intersections, angles) |
| `scipy.spatial.cKDTree` | Spatial nearest-neighbor lookups |
| `geopandas` | Spatial DataFrames for results |
