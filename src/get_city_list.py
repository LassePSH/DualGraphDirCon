import requests

# Set minimum population
min_pop = 2_000_000

# Query Overpass API
query = f"""
[out:json];
node["place"="city"]["population"];
out;
"""

response = requests.post(
    "https://overpass-api.de/api/interpreter",
    data=query
)

cities = response.json()['elements']

# Filter and print
for city in cities:
    pop_str = city['tags'].get('population', '0')
    # Remove spaces and non-digit characters
    pop_str = pop_str.replace(' ', '').replace(',', '')
    try:
        pop = int(pop_str)
        if pop >= min_pop:
            name = city['tags'].get('name:en')
            print(f"{name}: {pop:,}")
    except ValueError:
        continue