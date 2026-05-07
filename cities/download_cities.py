import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import osmnx as ox
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
from multiprocessing import Pool

CACHE_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'city_graphs')

COPENHAGEN_MUNICIPALITIES = [
    "Albertslund", "Ballerup", "Brøndby", "Dragør", "Frederiksberg",
    "Furesø", "Gentofte", "Gladsaxe", "Glostrup", "Greve",
    "Herlev", "Hvidovre", "Høje-Taastrup", "Ishøj", "Københavns",
    "Lyngby-Taarbæk", "Rudersdal", "Rødovre", "Tårnby", "Vallensbæk",
]


def _download_copenhagen():
    gdfs = []
    for name in COPENHAGEN_MUNICIPALITIES:
        try:
            gdfs.append(ox.geocode_to_gdf(f"{name} Kommune, Denmark"))
        except Exception as e:
            print(f'  Copenhagen: failed to fetch {name}: {e}')
    combined = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True))
    minx, miny, maxx, maxy = combined.geometry.union_all().bounds
    return ox.graph.graph_from_polygon(box(minx, miny, maxx, maxy))


def city_cache_path(city):
    safe = city.replace('/', '_').replace(' ', '_')
    return os.path.join(CACHE_DIR, safe + '.graphml')


def download_city(city):
    path = city_cache_path(city)
    if os.path.exists(path):
        print(f'{city} already cached, skipping')
        return
    try:
        if city.strip().lower().startswith('copenhagen'):
            G = _download_copenhagen()
        else:
            G = ox.graph.graph_from_place(city)
        ox.save_graphml(G, path)
        print(f'{city} downloaded and saved')
    except Exception as e:
        print(f'{city} failed: {e}')


def load_city(city):
    """Load a city graph from cache. Raises FileNotFoundError if not cached."""
    path = city_cache_path(city)
    if not os.path.exists(path):
        raise FileNotFoundError(f'No cached graph for {city!r}. Run download_cities.py first.')
    return ox.load_graphml(path)


if __name__ == '__main__':
    cities_file = os.path.join(os.path.dirname(__file__), 'cities.txt')
    with open(cities_file) as f:
        cities = [line.strip() for line in f if line.strip()]

    os.makedirs(CACHE_DIR, exist_ok=True)
    print(f'Downloading {len(cities)} cities to {os.path.abspath(CACHE_DIR)} ...')

    with Pool(90) as p:
        p.map(download_city, cities)

    print('Done.')
