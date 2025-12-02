from pyrosm import OSM, get_data
import os
import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
from shapely.geometry import Point, Polygon, MultiPoint, LineString, MultiLineString
from tqdm import tqdm
from shapely import geometry, ops
import shapely
import momepy
from my_tools import *
import shutup
shutup.please()
from multiprocessing import Pool

# OUTDATED

def get_degree_seq(city):
    fp = get_data(city)
    osm = OSM(fp)

    shape_df = osm.get_network(nodes=False)
    shape_df.crs = "epsg:4326"
    shape_df = shape_df.to_crs(3857)
    shape_df = shape_df.explode(index_parts=False)
    shape_df = gpd.GeoSeries(shape_df.unary_union.intersection(shape_df.unary_union), crs=shape_df.crs).explode()
    shape_df = gpd.GeoDataFrame(geometry=shape_df)
    G = momepy.gdf_to_nx(shape_df, approach="primal"); del shape_df
    G = clean_chains(G)
    points, lines = momepy.nx_to_gdf(G); del points
    G = momepy.gdf_to_nx(lines,approach='dual', multigraph=False,angles=False)
    G=new_angles(G)
    H, _, _ = merged_G_angle(G,tresh=10); del _
    degree_sequence = np.array([d for n, d in H.degree()]); del H
    return degree_sequence

def parallel(city):
    try:
        degree_sequence=get_degree_seq(city)
        np.savetxt(path+city+'.out', degree_sequence, delimiter=','); del degree_sequence
        print(city + ' succes!')
    except:
        print(city + ' failed')


if __name__ == '__main__':

    cities = ['Berlin', 'Aarhus', 'London', 'Barcelona', 'Boston', 'San Francisco', 'Tokyo', 'Paris', 'Rome', 
            'New York City', 'Sydney', 'Moscow', 'Amsterdam', 'Prague', 'Cairo', 'Bangkok', 'Dubai', 'Mumbai', 'Rio de Janeiro',
            'Toronto', 'Los Angeles', 'Singapore', 'Hong Kong', 'Seoul', 'Shanghai', 'Berlin', 'Cape Town', 'Buenos Aires',
            'Vienna', 'Florence', 'Dublin', 'Madrid', 'Athens', 'Stockholm', 'Dublin', 'Lisbon', 'Hanoi', 'Kuala Lumpur',
            'Oslo', 'Warsaw', 'Brussels', 'Budapest', 'Copenhagen', 'Edinburgh', 'Helsinki', 'Istanbul', 'Jerusalem',
            'Kyoto', 'Manila', 'Montreal', 'Nairobi', 'Ottawa', 'Panama City', 'Perth', 'Reykjavik', 'Santiago',
            'Vancouver', 'Wellington', 'Zurich', 'Abu Dhabi', 'Beijing', 'Brasília', 'Brisbane', 'Copenhagen', 'Doha',
            'Frankfurt', 'Geneva', 'Guangzhou', 'Hamburg', 'Jakarta', 'Kiev', 'Lima', 'Luxembourg City', 'Mexico City',
            'Monaco', 'Naples', 'Nice', 'Prague', 'Puerto Vallarta', 'San Diego', 'San Jose', 'San Juan', 'Santa Barbara',
            'São Paulo', 'Split', 'St. Petersburg', 'Tallinn', 'Valletta', 'Venice', 'Victoria', 'Vienna', 'Warsaw',
            'Washington D.C.', 'Wellington', 'Winnipeg', 'Yangon', 'Yerevan', 'Zagreb', 'Palermo', 'Palma de Mallorca',
            'Phnom Penh', 'Porto', 'Saint Petersburg', 'Salzburg', 'Sarajevo', 'Sofia', 'Tbilisi', 'Tirana', 'Ulaanbaatar',
            'Vilnius', 'Zanzibar City', 'Algiers', 'Antwerp', 'Asunción', 'Bangui', 'Belfast', 'Belgrade', 'Bishkek',
            'Bratislava', 'Bucharest', 'Cape Town', 'Caracas', 'Casablanca', 'Chisinau', 'Dakar', 'Dar es Salaam',
            'Dhaka', 'Douala', 'Freetown', 'Gaborone', 'Georgetown', 'Harare', 'Havana', 'Helsinki', 'Islamabad',
            'Kabul', 'Kampala', 'Kathmandu', 'Kigali', 'Kingston', 'Kolkata', 'Kuwait City', 'Libreville', 'Lilongwe',
            'Ljubljana', 'Lome', 'Luanda', 'Lusaka', 'Maseru', 'Maputo', 'Nassau', 'Niamey', 'Nouakchott',
            'Nukuʻalofa', 'Port-au-Prince', 'Port Louis', 'Port Moresby', 'Port of Spain', 'Port Vila', 'Praia', 'Pristina',
            'Quito', 'Rabat', 'Riga', 'San Salvador', 'Santo Domingo', 'Singapore', 'Suva', 'Tallinn', 'Tarawa', 'Tegucigalpa',
            'Thimphu', 'Tunis', 'Ulan Bator', 'Vaduz', 'Windhoek', 'Yamoussoukro','Copenhagen']



    path = '/work/lpsha/data/regime_v3/'

    existing = os.listdir(path)
    cities=[c.lower() for c in cities]
    existing=[ex[:-4] for ex in existing]

    cities_new = list(set(cities)-set(existing))    
    print(str(len(existing))  + ' out of '+str(len(cities)))

    with Pool(10) as p:
        p.map(parallel, cities_new)