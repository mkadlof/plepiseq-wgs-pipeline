#!/usr/bin/python

# TO DO

import  sys
from urllib.request import urlopen
from urllib.error import HTTPError
import urllib
import base64
import json
import numpy as np
import time
import os


API_TOKEN=open('/home/update/enterobase_api.txt').readlines()[0].rstrip()

def __create_request(request_str):
    base64string = base64.b64encode('{0}: '.format(API_TOKEN).encode('utf-8'))
    headers = {"Authorization": "Basic {0}".format(base64string.decode())}
    request = urllib.request.Request(request_str, None, headers)
    return request

DATABASE=sys.argv[1] # Nazwa bazy organizmu sentica, ecoli itd ..
CGNAME=sys.argv[2] # Nazwa ekwiwalentu bazy cg dla danego organizmu  w ecoli to cgMLST

# downloading entrire strain table
print('Downloading Strains table')
address = f'https://enterobase.warwick.ac.uk/api/v2.0/{DATABASE}/strains?my_strains=false&sortorder=asc&return_all=true&offset=0'
response = urlopen(__create_request(address))
data = json.load(response)
if os.path.exists('strains_table.npy'):
    # Wyciagamy informacje jakie barcode byly procesowane poprzednio
    old_strain_table = np.load('strains_table.npy', allow_pickle=True)
    old_strain_info = {}
    for strain in old_strain_table:
        old_strain_info[strain['strain_barcode']] = ''
    # Tutaj zapisujemy tylko strain_barcode ktore beda wykorzystane do przepytania baz STS i straindata
    strain_info= {}
    for strain in data['Strains']:
        if strain['strain_barcode'] not in old_strain_info.keys():
            strain_info[strain['strain_barcode']] = ''

else:
    # extract all strains from data first run
    strain_info= {}
    for strain in data['Strains']:
        strain_info[strain['strain_barcode']] = ''


# If straindata_table exists read it we will only update this dict with nove strains
if os.path.exists('straindata_table.npy'):
    straindata = np.load('straindata_table.npy', allow_pickle = True)
    straindata = straindata.item()
else:
    straindata = {}


# querying "strain" table in f'{step}'-sized chunks
# 150 seems to be optimal number , larger quieries are rejected
if len(strain_info.keys()) == 0:
    print('No new Strains in the Enterobase')
else:
    start = 0
    step = 160
    orignal_list = list(strain_info.keys())
    print(f'Downloading data for {len(orignal_list)} strains')
    if len(orignal_list) > step:
        for end in range(step,len(orignal_list), step):
            # create chunls
            print(f'Processing chunk: {start}:{end}')
            list_of_ids = orignal_list[start:end] # extract barcodes in chunks
            # build url link
            address2 = f"https://enterobase.warwick.ac.uk/api/v2.0/{DATABASE}/straindata?limit={step}&sortorder=asc&" + ('barcode={}&'*step).format(*list_of_ids) + "offset=0"
            response2 = urlopen(__create_request(address2))
            data2 = json.load(response2) # zwracamy tylko NAJMNIEJSZY ST z listyata2 = json.load(response2)
            # parse output basically we remove data from strain_info if the do not have required Sequence type
            straindata.update(data2['straindata'])
    
            start = end
        
            # We dont want to overburden the Enterobase API, downloading the database from scratch will take however few hours
            time.sleep(10)
    else:
        #original list is shorter than step size (we are probably doing an update and ww only need to dowload data for a few strains)
        end = 0

    # doliczmy brakujace id jesli nie pokryla ich petla
    if end != len(orignal_list):
        list_of_ids = orignal_list[start:len(orignal_list)]
        address2 = f"https://enterobase.warwick.ac.uk/api/v2.0/{DATABASE}/straindata?limit={step}&sortorder=asc&" + ('barcode={}&'*len(list_of_ids)).format(*list_of_ids) + "offset=0"
        response2 = urlopen(__create_request(address2))
        data2 = json.load(response2) # zwracamy tylko NAJMNIEJSZY ST z listyata2 = json.load(response2)
        # parse output basically we remove data from strain_info if the do not have required Sequence type
        straindata.update(data2['straindata'])


# Extracting known STs
print('Analyzing STs')
sts_data = {}
known_STs = []
old_known_STs = []

if os.path.exists('sts_table.npy'):
    # wczytyjemy stare dane do sts_data
    sts_data = np.load('sts_table.npy', allow_pickle=True)
    sts_data = sts_data.item()
    #  w tabeli STs wartosc ST jest str ale w tabeli straindata jest int-em
    # jako ze bedzie powrownywac sie z straindata wracamy do typu int tymczasowo
    #old_known_STs = list(map(int, list(sts_data.keys())))


# iterujemy po straindata
# ktore na tym etapie zawiera pelna baze entero (zarowno w trybie update jak i de novo)
for klucz,wartosc in straindata.items():
    try:
        for scheme in wartosc['sts']:
            if CGNAME in scheme.values():
                if scheme['st_id'] >= 1 and str(scheme['st_id']) not in sts_data.keys():
                    # w tabeli straindata st_id jest intem
                    known_STs.append(scheme['st_id'])
                else:
                    pass
                    # CZasami w bazie sa ST niejsze od 0
                    # patrz przyklad SAL_AC3218AA
    except:
        # Stare strainy 
        # nie zawsze maja pole 'sts'
        pass
#known STs zawiera tylko informacje o STs do pobrania

known_STs = sorted(list(set(known_STs)))

# enterobase sugeruje sciagac max 500 wartosci, ograniczmy sie do 400 + sleep 10s

if len(known_STs) == 0:
    print('No new pHierCC data to download')
else:
    print(f'Need to download {len(known_STs)} novel STs')
    start = 0
    step = 400
    if len(known_STs) > step:
        for end in range(step,len(known_STs), step):
            # create chunls
            print(f'Processing chunk: {start}:{end}')
            list_of_ids = known_STs[start:end]
            address3 = f'https://enterobase.warwick.ac.uk/api/v2.0/{DATABASE}/{CGNAME}/sts?limit=500&scheme={CGNAME}&st_id=' + f'{list_of_ids[0]}' + ('%2C{}'*(len(list_of_ids)-1)).format(*list_of_ids[1:]) + '&offset=0'
            response3 = urlopen(__create_request(address3))
            data3 = json.load(response3)
            for wpis in data3['STs']:
                sts_data[wpis['ST_id']] = wpis['info']['hierCC']
            start = end
            time.sleep(10)
    else:
        end  = 1

    if end != len(known_STs):
        list_of_ids = known_STs[start:len(known_STs)]
        address3 = f'https://enterobase.warwick.ac.uk/api/v2.0/{DATABASE}/{CGNAME}/sts?limit=500&scheme={CGNAME}&st_id=' + f'{list_of_ids[0]}' + ('%2C{}'*(len(list_of_ids)-1)).format(*list_of_ids[1:]) + '&offset=0'
        response3 = urlopen(__create_request(address3))
        data3 = json.load(response3)
        for wpis in data3['STs']:
            sts_data[wpis['ST_id']] = wpis['info']['hierCC']

# Save all the tables
np.save('strains_table.npy', data['Strains'], allow_pickle=True, fix_imports=True)
np.save('straindata_table.npy', straindata, allow_pickle=True, fix_imports=True)
np.save('sts_table.npy', sts_data, allow_pickle=True, fix_imports=True)
