import sys
sys.path.append('/home/lpsha/s154446/fractality/dual_graph_algo/')
import dual_conti
import numpy as np
from multiprocessing import Pool
import osmnx as ox


def get_degree_seq(city):
    G = ox.graph.graph_from_place(city)
    gdf_merged, _, _, _ = dual_conti.get_dual_dir_con(t_buffer=10, a_threshold=20, data=G, enforce_degree2=False)
    degree_sequence = np.array(gdf_merged.degree)
    return degree_sequence

def parallel(city):
    try:
        degree_sequence=get_degree_seq(city)
        np.savetxt(path+city+'.out', degree_sequence, delimiter=','); del degree_sequence
        print(city + ' succes!')
    except:
        print(city + ' failed')


if __name__ == '__main__':

    cities = ["London, England, UK", "Copenhagen, Denmark", "Berlin, Germany", "Aarhus, Denmark", "Barcelona, Spain", "Boston, Massachusetts, USA",
             "San Francisco, California, USA", "Tokyo, Japan", "Paris, France", "Rome, Italy", 
            "New York City, New York, USA", "Sydney, Australia", "Moscow, Russia", "Amsterdam, Netherlands", "Prague, Czechia", "Cairo, Egypt", "Bangkok, Thailand", "Dubai, United Arab Emirates", "Mumbai, India", "Rio de Janeiro, Brazil",
            "Toronto, Canada", "Los Angeles, California, USA", "Singapore", "Hong Kong, China", "Seoul, South Korea", "Shanghai, China", "Cape Town, South Africa", "Buenos Aires, Argentina",
            "Vienna, Austria", "Florence, Italy", "Dublin, Ireland", "Madrid, Spain", "Athens, Greece", "Stockholm, Sweden", "Lisbon, Portugal", "Hanoi, Vietnam", "Kuala Lumpur, Malaysia",
            "Oslo, Norway", "Warsaw, Poland", "Brussels, Belgium", "Budapest, Hungary", "Edinburgh, Scotland, UK", "Helsinki, Finland", "Istanbul, Turkey", 
            "Kyoto, Japan", "Manila, Philippines", "Montreal, Canada", "Nairobi, Kenya", "Ottawa, Canada", "Panama City, Panama", "Perth, Australia", "Reykjavik, Iceland", "Santiago, Chile",
            "Vancouver, Canada", "Wellington, New Zealand", "Zurich, Switzerland", "Abu Dhabi, United Arab Emirates", "Beijing, China", "Brasília, Brazil", "Brisbane, Australia", "Doha, Qatar",
            "Frankfurt, Germany", "Geneva, Switzerland", "Guangzhou, China", "Hamburg, Germany", "Jakarta, Indonesia", "Kiev, Ukraine", "Lima, Peru", "Luxembourg City, Luxembourg", "Mexico City, Mexico",
            "Monaco", "Naples, Italy", "Nice, France", "Puerto Vallarta, Mexico", "San Diego, California, USA", "San Jose, California, USA", "San Juan, Puerto Rico", "Santa Barbara, California, USA",
            "São Paulo, Brazil", "Split, Croatia", "St. Petersburg, Russia", "Tallinn, Estonia", "Valletta, Malta", "Venice, Italy", "Victoria, Canada", 
            "Washington, D.C., USA", "Winnipeg, Canada", "Yangon, Myanmar", "Yerevan, Armenia", "Zagreb, Croatia", "Palermo, Italy", "Palma de Mallorca, Spain",
            "Phnom Penh, Cambodia", "Porto, Portugal", "Salzburg, Austria", "Sarajevo, Bosnia and Herzegovina", "Sofia, Bulgaria", "Tbilisi, Georgia", "Tirana, Albania", "Ulaanbaatar, Mongolia"]

    print("Number of cities:", len(cities))

    path = '/home/lpsha/s154446/fractality/data/city_degrees/t10_a20/' # out path

    with Pool(20) as p:
        p.map(parallel, cities)