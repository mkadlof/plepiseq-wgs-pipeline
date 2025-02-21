import pandas as pd
import geodatasets
import geopandas
import folium
import datetime
import calendar
import matplotlib.colors as mcolors
from folium.plugins import TimeSliderChoropleth
import numpy as np
import sys
import itertools as it
import time

# by deafult data are split into 20 bins, each bi with its own color
# to change that replace n_col value with smt else

# by default the gradient is from light to dark blue
# to change that replace color codes in light_blue and dark blue variables

# The script requires only two parameters input (enterobase historical data obtained within the pipeline
# with "extract_historical_data_enterobase" module
# and the prefix for the name of the output file

dataframe = pd.read_csv(sys.argv[1], sep = "\t", dtype = 'str', na_values = 'None')
wynik = dataframe[['Country', 'Year']].groupby(['Country', 'Year']).size()
max_value = wynik.max()

max_year = dataframe.Year.dropna().astype(int).max() + 1
min_year = dataframe.Year.dropna().astype(int).min() - 1

white= "#ffffff"
light_blue = "#dcf3fa"
dark_blue = "#164452"
cmap = mcolors.LinearSegmentedColormap.from_list("blue_gradient", [light_blue, dark_blue])

if dataframe.Country.size < 1:
    # empty data but we need to produce something for nexflow
    n_col = 2
elif dataframe.Country.size <= 20:
    n_col = dataframe.Country.size + 1
else:
    n_col = 20
gradient_colors = [mcolors.to_hex(cmap(i / (n_col - 1))) for i in range(n_col)]

# geojeson file within container
dane = geopandas.read_file('/data/ne_10m_admin_0_countries.geojson')

style_dict = {}

for country in dataframe.Country.dropna().unique():
    for year in range(min_year,max_year, 1):
        try:
            indeks = dane[(dane.SOVEREIGNT == country) & (dane.scalerank == 0)].index[0]
            if str(indeks) not in style_dict.keys():
                style_dict[str(indeks)] = {}
            czas = datetime.datetime(int(year),1, 1, 0, 0, 0)
            czas = calendar.timegm(czas.timetuple())
            if str(czas) not in  style_dict[str(indeks)].keys():
                style_dict[str(indeks)][str(czas)]= {}
            try: 
                ilosc=wynik[(str(country), str(year))]
                color = gradient_colors[np.where(np.histogram(wynik[(str(country), str(year))], bins=n_col, range=[0, max_value])[0] == 1)[0][0]]
                opacity = 1
            except KeyError:
                #print(f'Brak wpisu dla  {country} {year}')
                color = "#ffffff"
                opacity = 0
            style_dict[str(indeks)][str(czas)] = {'color' :color, 'opacity':opacity}

        except:
            print(f'Error for {country} {year}')

if len(style_dict.keys()) == 0:
    # Nie ma danych w enterobase o tym cluster, ale porbka jest z polski
    style_dict['76'] = {}
    style_dict['76'][str(int(time.time()))] = {'color' : light_blue, 'opacity':1}

m = folium.Map([48, 2], zoom_start=5)
TimeSliderChoropleth(dane.to_json(), styledict=style_dict).add_to(m)
m.save(f'{sys.argv[2]}.html')


