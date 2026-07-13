# Extend the city list for wider street-type diversity

**Date:** 2026-07-13
**Status:** approved

## Problem

The analysis set (`cities/cities.txt`, 96 cities; `data/city_external/city_info.csv`,
97 rows) is morphologically lopsided. Current `street_pattern` distribution:

| pattern | count |
|---------|-------|
| organic | 58 |
| grid    | 25 |
| mixed   | 5 |
| radial  | 5 |
| planned | 4 |

Geographically it is Europe-heavy (46/97) with thin Oceania (2) and Africa (6).
The user wants "a wider range of types of streets" — richer coverage of the
under-represented pattern buckets *and* regional traditions, along four axes:
under-represented patterns, informal/unplanned morphologies, extreme/distinctive
grids, and geographic breadth.

## Approach

Add **40 curated cities** (mixed sizes, per user) spanning the four axes, none
duplicating the existing 96. Each city is appended to `cities.txt` and gets a
full 14-column row in `city_info.csv`. All 40 names were validated to resolve to
an OSM `(Multi)Polygon` boundary via `osmnx.geocode_to_gdf` before selection.

### The 40 cities

**Radial / planned / modernist (8)**
Canberra AU · Brasília BR · New Delhi IN · Karlsruhe DE · Adelaide AU ·
Chandigarh IN · La Plata AR · Palmanova IT

**Informal / medina / non-European organic (10)**
Marrakesh MA · Fez MA · Tunis TN · Isfahan IR · Kathmandu NP · Medellín CO ·
Mumbai IN · Kampala UG · Lahore PK · Kigali RW

**Extreme / distinctive grids (8)**
Chicago US · Manhattan (New York) US · Salt Lake City US · Melbourne AU ·
Turin IT · Savannah US · Beijing CN · Johannesburg ZA

**Geographic breadth (14)**
Dakar SN · Dar es Salaam TZ · Cape Town ZA · Tashkent UZ · Almaty KZ ·
Samarkand UZ · Colombo LK · Auckland NZ · Quito EC · Guatemala City GT ·
Tehran IR · Chengdu CN · Istanbul TR · Kolkata IN

### Resulting balance

- `street_pattern`: mixed 5→17, radial 5→10, planned 4→6, grid 25→36, organic 58→68.
- Continents: Oceania 2→6, Africa 6→15, plus Central/South Asia and Latin America.
- Grid end stressed with hyper-regular (Chicago/Manhattan/SLC), Roman-castrum
  (Turin), ward-square (Savannah), and modernist-sector (Chandigarh).

## Name resolution fix

Five names (Canberra, Mumbai, Lahore, Johannesburg, Cape Town) resolve to a
polygon only with `which_result=1`; with osmnx's default (`which_result=None`,
`limit=50`) Nominatim returns a node first and `graph_from_place` raises
"did not geocode … to a (Multi)Polygon". `cities/download_cities.py` currently
calls `graph_from_place(city)` with the default, so those five would fail to
download.

Fix: in `download_city`, on failure retry `graph_from_place(city, which_result=1)`
before giving up. Minimal, matches the existing try/except-and-print pattern, and
makes any future quirky name robust.

## city_info.csv rows

All 14 columns are populated per city (no blanks), using compiled ~2020–2023
estimates consistent with the existing rows and the README's stated caveats
(city-proper vs metro boundaries vary; GDP is country-level; `street_pattern` is
a dominant-type simplification). `population_density_per_km2` is derived as
`round(population / area_km2)`.

## Deliverables

1. Append 40 lines to `cities/cities.txt`.
2. Append 40 rows to `data/city_external/city_info.csv`.
3. Add the `which_result=1` fallback to `cities/download_cities.py`.
4. Update `data/city_external/README.md` count 97 → 137.
5. Commit.

## Out of scope

Downloading the graphs / running DGDC on the new cities. `download_cities.py`,
`regime.py`, etc. pick them up from `cities.txt` on the next run; no code beyond
the fetch fallback changes.
