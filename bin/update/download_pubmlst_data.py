import requests

import subprocess
import sys
import time
import numpy as np
import os
from multiprocessing import Pool


def execute_command(polecenie: str):
    """

    :param polecenie: str, polecenie do wykonania w powloce
    :return: Polecenie ma zwrocic zwykle jakis plik, wiec sama funkcja nie zwraca nic
    """
    try:
        wykonanie = subprocess.Popen(polecenie, shell=True, stdout=subprocess.PIPE)
        wykonanie.communicate()
        return True
    except:
        return False

def download_pubmlst_profiles(profile_links, start, end, worker_id = -1):
    i = 0
    local_dict = {}
    for profile in profile_links[start:end]:
        profile_data = requests.get(profile).json()
        local_dict[str(profile_data['cgST'])] = {}
        for level, level_data in profile_data['classification_schemes'].items():
            level_number = level.split('_')[-1]
            level_str = f'd{level_number}'
            local_dict[str(profile_data['cgST'])][level_str] = level_data['group']['group']
        if i % 1000 == 0:
            print(f"Analyzing profile {profile_data['cgST']} in worker {worker_id} " 
                  f"on {time.strftime('%Y-%m-%d %H:%M', time.gmtime())}")

        i += 1
    return local_dict
def download_pubmlst_entries(list_of_isolates, start, end, DATABASE,  DATABASE2, worker_id=-1):
    missing_link = open(f"missing_isolates_{worker_id}.txt", "w", buffering=1)
    local_dict = {}
    i = 0
    for isolate in list_of_isolates[start:end]:
        if i % 1000 == 0:
            print(f"Processing link no. {i} in worker {worker_id} "
                  f"on {time.strftime('%Y-%m-%d %H:%M', time.gmtime())}")
            time.sleep(4)
        i += 1
        try:
            isolate_info = requests.get(isolate).json()
            isolate_key = str(isolate_info['provenance']['id'])

            # Basic info should always be available
            isolate_id = isolate_info['provenance']['isolate']
            isolate_country = isolate_info['provenance']['country']
            isolate_date_entered = isolate_info['provenance']['date_entered']
            if 'year' in isolate_info['provenance'].keys():
                isolate_date_year = isolate_info['provenance']['year']
            else:
                isolate_date_year = isolate_info['provenance']['date_entered'].split('-')[0]

            # fill dict
            local_dict[isolate_key] = {}
            local_dict[isolate_key]['isolate_id'] = isolate_id
            local_dict[isolate_key]['country'] = isolate_country
            local_dict[isolate_key]['date_entered'] = isolate_date_entered
            local_dict[isolate_key]['year'] = isolate_date_year

            # check if there is sequencing info
            try:
                local_dict[isolate_key]['sequencing'] = isolate_info['isolate_info']['biosample_accession'][1]
            except KeyError:
                # This data is only available for some samples
                local_dict[isolate_key]['sequencing'] = 'No data'

            # add data regarding schemes
            local_dict[isolate_key]['sts'] = []  # STs are stored as list of dicts
            # isolates_dict[isolate_key]['hiercc'] = {} There is no nned to declare it now

            MLST_data = 0
            cgMLST_data = 0
            hierCC_data = 0
            if 'schemes' in isolate_info.keys():
                # some isolate have no 'schemes' entry
                for scheme in isolate_info['schemes']:
                    #  If we can assign MLST to an isolate MLST_data is set to 1
                    #  otherwise below we assign it a dummy 'unk' value
                    if scheme['description'] == DATABASE2:
                        if 'fields' in scheme.keys():
                            # isolate can have a scheme but without assigne ST/cgST value
                            if 'ST' in scheme['fields'].keys():
                                MLST_data = 1
                                local_dict[isolate_key]['sts'].append(
                                    {'scheme_name': DATABASE2, 'st_id': scheme['fields']['ST']})

                    elif scheme['description'] == DATABASE:
                        level_0 = 'unk'
                        if 'fields' in scheme.keys():
                            # some isolate have partial cgST assignment
                            # for example the entry is present but without cgST itself (like id 78351)
                            # in some cases cgST is not as int but a list cgST = [number1, number2 ...] and so on

                            if 'cgST' in scheme['fields'].keys():

                                if isinstance(scheme['fields']['cgST'], int) or isinstance(scheme['fields']['cgST'],
                                                                                           str):
                                    cgMLST_data = 1
                                    local_dict[isolate_key]['sts'].append(
                                        {'scheme_name': DATABASE, 'st_id': scheme['fields']['cgST']})
                                    level_0 = scheme['fields']['cgST']
                                elif isinstance(scheme['fields']['cgST'], list):
                                    cgMLST_data = 1
                                    local_dict[isolate_key]['sts'].append(
                                        {'scheme_name': DATABASE, 'st_id': scheme['fields']['cgST'][0]})
                                    level_0 = scheme['fields']['cgST'][0]

                        if 'classification_schemes' in scheme.keys():
                            # classification_schemes can only present for schemes related to cgMLSt
                            # however clustering is independent from presence/abscence of cgST assignment itself
                            # so we can know to which cluster and isolate belongs, but not what its cgST is ...
                            # thus the 'd0' key might be unk, yet other keys are ok
                            level_5 = 'unk'
                            level_10 = 'unk'
                            level_25 = 'unk'
                            level_50 = 'unk'
                            level_100 = 'unk'
                            level_200 = 'unk'

                            try:
                                # just in case some keys are missing from pubmlst
                                if 'Cjc_cgc2_5' in scheme['classification_schemes'].keys():
                                    level_5 = scheme['classification_schemes']['Cjc_cgc2_5']['groups'][0]['group']
                                if 'Cjc_cgc2_10' in scheme['classification_schemes'].keys():
                                    level_10 = scheme['classification_schemes']['Cjc_cgc2_10']['groups'][0]['group']
                                if 'Cjc_cgc2_25' in scheme['classification_schemes'].keys():
                                    level_25 = scheme['classification_schemes']['Cjc_cgc2_25']['groups'][0]['group']
                                if 'Cjc_cgc2_50' in scheme['classification_schemes'].keys():
                                    level_50 = scheme['classification_schemes']['Cjc_cgc2_50']['groups'][0]['group']
                                if 'Cjc_cgc2_100' in scheme['classification_schemes'].keys():
                                    level_100 = scheme['classification_schemes']['Cjc_cgc2_100']['groups'][0]['group']
                                if 'Cjc_cgc2_200' in scheme['classification_schemes'].keys():
                                    level_200 = scheme['classification_schemes']['Cjc_cgc2_200']['groups'][0]['group']
                                local_dict[isolate_key]['hiercc'] = {'d0': level_0, 'd5': level_5, 'd10': level_10,
                                                                     'd25': level_25, 'd50': level_50,
                                                                     'd100': level_100, 'd200': level_200}
                                hierCC_data = 1
                            except:
                                # hierCC_data will be 0 and below we put "unk" as values
                                pass
            else:
                pass
                # we put dummy values for 'sts' and 'hiercc'
                # print(f'No schemes for id {i}')

            # if a specific type of data is missing add it with dummy values
            if MLST_data == 0:
                # put dummy values for that isolate
                local_dict[isolate_key]['sts'].append({'scheme_name': DATABASE2, 'st_id': 'unk'})

            if cgMLST_data == 0:
                local_dict[isolate_key]['sts'].append({'scheme_name': DATABASE, 'st_id': 'unk'})

            if hierCC_data == 0:
                local_dict[isolate_key]['hiercc'] = {'d0': 'unk', 'd5': 'unk', 'd10': 'unk', 'd25': 'unk',
                                                     'd50': 'unk', 'd100': 'unk', 'd200': 'unk'}

        except:
            missing_link.write(f"{isolate}\n")
            print(f'Error when downloading data for isolate {isolate} for worker {worker_id}')
    missing_link.close()
    return local_dict


# Data for hiercc 
SPEC = "pubmlst_campylobacter_seqdef"
DATABASE = 'C. jejuni / C. coli cgMLST v2'

# Data for straindata
SPEC2 = "pubmlst_campylobacter_isolates"  # isolates database
DATABASE2 = 'MLST'  # additional database for isolates

CURRENT_DATE=time.strftime("%Y-%m-%d", time.gmtime())

# Set the maximum number of cpus to 4 to avoid to many simultaneous requests
cpus = int(sys.argv[1])
if cpus > 4:
    cpus = 4

#  we need to  find schema link (there are many schemas like MLST , cgMLST, rMLST ...)
#  cgMLST have assosiciated scheme_classification data (phiercc equivalent from enterobase)
scheme_link = ''
scheme_table = requests.get('https://rest.pubmlst.org/db/' f'{SPEC}' + '/schemes')
for scheme in scheme_table.json()['schemes']:
    if DATABASE == scheme['description']:
        scheme_link = scheme['scheme']

# download 
# check if this is an update or set up new dict
if os.path.exists('sts_table.npy'):
    hiercc_dict = np.load('sts_table.npy', allow_pickle=True)
    hiercc_dict = hiercc_dict.item()
else:
    hiercc_dict = {}


if os.path.exists('timestamp'):
    previous_update = open('timestamp').readlines()[0]
else:
    previous_update = '1990-01-01'


profiles = requests.get(scheme_link + '/profiles', {'updated_after': previous_update})
record_list = profiles.json()['records']
print(f'Downloading {record_list} profiles available since {previous_update}')

profiles = requests.get(scheme_link + '/profiles', {'page_size': record_list}).json()

pool = Pool(cpus)
lista_indeksow = []
start = 0
step = len(profiles['profiles']) // cpus

for i in range(cpus):
    end = start + step
    if i == (cpus - 1) or end > len(profiles['profiles']):
        # if this is a last cpu, or one has passed the last element of a list, use length.
        end = len(profiles['profiles'])
    lista_indeksow.append([start, end])
    start = end

# if there are less than 20 profiles, we use a single cpu
if len(profiles['profiles']) < 20:
    lista_indeksow = [[0, len(profiles['profiles'])]]

jobs = []
worker_id = 1  # we use that only for logging
for start, end in lista_indeksow:
    jobs.append(pool.apply_async(download_pubmlst_profiles, (profiles['profiles'], start, end, worker_id)))
    worker_id += 1
pool.close()

# update dict with workers results
for job in jobs:
    hiercc_dict = {**hiercc_dict, **job.get()}
pool.join()

np.save('sts_table.npy', hiercc_dict, allow_pickle=True, fix_imports=True)
# Get/Update strain data

# Repeat the same code for straindata
if os.path.exists('straindata_table.npy'):
    isolates_dict = np.load('straindata_table.npy', allow_pickle=True).item()
else:
    isolates_dict = {}

all_isolates = requests.get(f'https://rest.pubmlst.org/db/{SPEC2}/isolates',
                            params={'updated_after':previous_update})
all_isolates_number = all_isolates.json()['records']
all_isolates = requests.get(f'https://rest.pubmlst.org/db/{SPEC2}/isolates',
                            params={'page_size': all_isolates_number, 'updated_after':previous_update}).json()

print(f'Downloading data for {all_isolates_number} isolates available since {previous_update}')

pool = Pool(cpus)
lista_indeksow = []
start = 0
step = len(all_isolates['isolates']) // cpus

for i in range(cpus):
    end = start + step
    if i == (cpus - 1) or end > len(all_isolates['isolates']):
        # if this is a last cpu or yoy passed the last element of a list, use as end len(lista_alleli)
        end = len(all_isolates['isolates'])
    lista_indeksow.append([start, end])
    start = end

# if there are less than 20 isolates, just use a single cpu
if len(all_isolates['isolates']) < 20:
    lista_indeksow = [[0, len(all_isolates['isolates'])]]

jobs = []
worker_id = 1
for start, end in lista_indeksow:
    jobs.append(pool.apply_async(download_pubmlst_entries, (all_isolates['isolates'], start, end, DATABASE, DATABASE2, worker_id)))
    worker_id += 1
pool.close()

# update dict with workers results
for job in jobs:
    isolates_dict = {**isolates_dict, **job.get()}
pool.join()

#  Check all the files named missing_isolates_{worker_id}.txt for isolates that were skipped
#  and try to download data again
time.sleep(60)
lista_plikow = [x for x in os.listdir(".") if "missing_isolates" in x]
lista_linkow = []

for plik in lista_plikow:
    lista_linkow.extend([x.rstrip() for x in open(plik).readlines()])
    os.system(f"rm {plik}")

print(f"Second attempt to download {len(lista_linkow)} entries")
if len(lista_linkow) > 0:
    additional_data = download_pubmlst_entries(list_of_isolates = lista_linkow,
                                               start=0,
                                               end=len(lista_linkow),
                                               DATABASE=DATABASE,
                                               DATABASE2=DATABASE2,
                                               worker_id=-1)
    isolates_dict = {**isolates_dict, **additional_data}

np.save('straindata_table.npy', isolates_dict, allow_pickle=True, fix_imports=True)

# set new last update date
with open('timestamp', 'w') as f:
    f.write(CURRENT_DATE)