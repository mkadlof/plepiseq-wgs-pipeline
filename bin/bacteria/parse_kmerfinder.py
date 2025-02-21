# micro script to parse kmerfinder output stored as a jason
import json
import sys

# sys,argv[1] is a json file
# sys.aegv[2] is either "genus" or "species"

full_data = json.load(open(sys.argv[1]))
if sys.argv[2] == 'genus':
    for klucz, wartosc in full_data['kmerfinder']['results']['species_hits'].items():
        print(wartosc['Taxonomy'].split(';')[-1].lstrip().split(" ")[0])
        break # just in case there are more than 2 species
elif sys.argv[2] == 'species':
    for klucz, wartosc in full_data['kmerfinder']['results']['species_hits'].items():
        print(wartosc['Taxonomy'].split(';')[-1].lstrip())
        break
else:
    pass
