import pandas as pd
import geopandas
import numpy as np
import sys
import itertools as it
import time
import plotly.express as px


min_year_to_plot = int(sys.argv[3]) # only plot data fter this year  

dataframe = pd.read_csv(sys.argv[1], sep = "\t", dtype = 'str', na_values = 'None')
wynik = dataframe[['Country', 'Year']].groupby(['Country', 'Year']).size()
max_value = wynik.max()

if dataframe.dropna().Year.astype(int).size == 0:
    min_year = min_year_to_plot
    max_year = min_year_to_plot + 1
else:
    max_year = dataframe.Year.dropna().astype(int).max() + 1
    min_year = dataframe.Year.dropna().astype(int).min() - 1

    if min_year < min_year_to_plot:
        min_year = min_year_to_plot

    if max_year < min_year_to_plot:
        max_year = min_year_to_plot + 1

# geojeson file within container
dane = geopandas.read_file('/data/ne_10m_admin_0_countries.geojson')
dane.geometry = dane.geometry.simplify(0.1) # reduces complexity of polygons, thus borders might look a bit off on larger zoom levels
# but it makes html MUCH smaller. The smaller the number the more accurate are data for borders

indexes = [] # index of a count in geojeson
years = []
val = [] # number of cases


year_order = []
for country in dataframe.Country.dropna().unique():
    for year in range(min_year,max_year, 1):
        year_order.append(year)
        try: 
            indeks = dane[(dane.SOVEREIGNT == country) & (dane.scalerank == 0)].index[0]
            ilosc=wynik[(str(country), str(year))]
            indexes.append(str(indeks))
            years.append(year)
            val.append(ilosc)
        except:
            # if a given countr in a given year has no cases we will ommit it
            # plotly does not like it but we ensure correct data visualization by creating year_order list

            # in case we keep 0 a larger number of polygons will be drawn on a map
            # and will SIGNIFICANTLY increase html size
            pass


if len(indexes) == 0:
    # no data in entorabase but we nned to draw something
    year_order = [2023, 2024]
    style_df = pd.DataFrame({'indexes':['76', '76'], 'year':[2023, 2024], 'wartosc':[0, 1]})
    max_val = 2
else:
    style_df = pd.DataFrame({'indexes':indexes, 'year':years, 'wartosc':val})
    if style_df.wartosc.max() > 100:
        max_val = 100
    else:
        max_val = style_df.wartosc.max()+1

# keep only indexes and polygons from geojeson
tmp_geojeson = dane.geometry

# do not show data prior to 2009
style_df = style_df[style_df.year >= 2009]

fig = px.choropleth_mapbox(style_df, \
        geojson=tmp_geojeson, \
        color='wartosc', \
        locations = 'indexes', \
        animation_frame="year", \
        range_color=(0, max_val),   \
        color_continuous_scale='Blues', \
        zoom=3, \
        center = {"lat": 48.0902, "lon": 2.7129}, \
        mapbox_style="open-street-map", \
        opacity = 0.5, \
        category_orders = {'year': year_order}, \
        labels = {'wartosc':'ilosc przypadkow'})

# save the data
with open(f'{sys.argv[2]}.html', 'w') as f:
    f.write(fig.to_html(full_html=False, include_plotlyjs=True))

