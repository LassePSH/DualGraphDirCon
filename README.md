# fractality

Street network morphology research project studying urban street patterns through dual graph continuity analysis and generative space-filling models.

## Overview

This project provides two main components:

1. **Dual Graph Directional Continuity (DGDC)** — analyzes real-world street networks by grouping street segments into continuous "axes" based on angular continuity
2. **Space-Filling Street Model** — generatively produces fractal-like street networks through iterative subdivision

## Installation

```bash
pip install -e .
```

Requires Python ≥ 3.9.

## Usage

### Dual Graph Analysis

```python
import osmnx as ox
from dgdc import get_dual_dir_con

G = ox.graph_from_place("Amsterdam, Netherlands", network_type="drive")

gdf_merged, H, shape_exploded_df, lines = get_dual_dir_con(
    t_buffer=10,       # buffer radius (m) around intersection touch points
    a_threshold=20,    # max angle (degrees) to consider streets continuous
    data=G,
    enforce_degree2=False
)
```

**Pipeline:**
1. Fetch OSM street network via `osmnx`
2. Merge degree-2 nodes into single edges (`clean_chains`)
3. Compute directional angles at intersections using touch-point buffers (`new_angles`)
4. Merge dual-graph nodes with similar angles below `a_threshold` (`merged_G_angle`)
5. Optionally enforce max degree 2 (`split_until_degree_2`)

### Space-Filling Street Model

```python
import numpy as np
from space_filling_model.space_filling import run_model, get_geodataframe

# Model 1 — basic midpoint subdivision
gdf = get_geodataframe(model=1, seed=42)

# Model 2 — biased subdivision with angle and area constraints
segs = run_model(model=2, seed=42, min_angle=np.pi/4)
```

**Models:**
- **Model 1 (basic):** Select segment proportional to length, split at midpoint, connect to nearest visible intersection
- **Model 2 (biased):** Same as Model 1 with constraints: polygon area `A(r) > 0.05 * exp(-1/r)` and all intersection angles > π/4

### City-Scale Analysis

Compute dual graph degree sequences for ~84 cities in parallel:

```bash
python cities/regime.py
```

Results are saved to `data/city_degrees/t10_a20/<city>.out`.

## Project Structure

```
dgdc/               # Core DGDC package
  dual_conti.py     # Main algorithm
  __init__.py
space_filling_model/
  space_filling.py  # Generative street model
cities/
  regime.py         # Batch city analysis (84 cities, 20 parallel processes)
  cities.txt        # List of cities
notebooks_lasse/    # Exploratory Jupyter notebooks
data/               # Output data
```

## Dependencies

| Package | Purpose |
|---------|---------|
| `osmnx` | Fetch OSM street networks |
| `momepy` | Urban morphology (roundabout simplification, primal/dual graph conversion) |
| `networkx` | Graph construction and analysis |
| `shapely` | Geometric operations |
| `scipy` | Spatial nearest-neighbor lookups (`cKDTree`) |
| `geopandas` | Spatial DataFrames |
| `numpy` | Numerical operations |
