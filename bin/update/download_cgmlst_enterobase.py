from urllib.request import urlopen
from urllib.error import HTTPError
import urllib
import base64
import json
import os
import subprocess
import time
from Bio import SeqIO
import sys
from multiprocessing import Pool

API_TOKEN = open('/home/update/enterobase_api.txt').readlines()[0].rstrip()
DATABASE = sys.argv[1]  # ecoli, yersinia, mcatarrhalis API 2.0
scheme_name = sys.argv[2]  # according to API 2.0
scheme_dir = sys.argv[3]  # take a look inside https://enterobase.warwick.ac.uk//schemes/ for names of specific schemes
cpus = int(sys.argv[4])  # Liczba CPUs uzywanych prz tworzeniu bazy dla blasta


def __create_request(request_str):
    base64string = base64.b64encode('{0}: '.format(API_TOKEN).encode('utf-8'))
    headers = {"Authorization": "Basic {0}".format(base64string.decode())}
    request = urllib.request.Request(request_str, None, headers)
    return request


def execute_command(polecenie: str):
    """

    :param polecenie: str, polecenie do wykonania w powloce
    :return: Polecenie ma zwrocic zwykle jakis plik, wiec sama funkcja nie zwraca nic
    """
    wykonanie = subprocess.Popen(polecenie, shell=True, stdout=subprocess.PIPE)
    wykonanie.communicate()
    return True


def run_blast(lista_loci, start, end):
    # index the file
    for locus in lista_loci[start:end]:
        execute_command(f"gunzip {locus}.fasta.gz")
        execute_command(f"makeblastdb -in {locus}.fasta -dbtype nucl")
    return True

def get_feader(file_path):
    """
    Simple function to extract allele number from headears int the fasta file.
    :return: dict, only one key - allele name, value list of ints (allel versions)
    """
    dictionary = {}
    allele_name=os.path.basename(file_path).replace('.fasta','')
    dictionary[allele_name] = []
    for seq_record in SeqIO.parse(file_path, "fasta"):
        dictionary[allele_name].append(seq_record.id.rsplit('_')[1])

    return dictionary

address = f'https://enterobase.warwick.ac.uk/api/v2.0/{DATABASE}/{scheme_name}/loci?limit=10000&scheme={scheme_name}&offset=0'
lista_loci = []
try:
    response = urlopen(__create_request(address))
    data = json.load(response)
    i = 0
    for locus in data['loci']:
        if i % 100 == 0:
            time.sleep(1)
        #print(locus['locus'])
        i += 1
        locus_link = f"https://enterobase.warwick.ac.uk//schemes/{scheme_dir}/{locus['locus']}.fasta.gz"
            # download the locus_address
        response_locus = urlopen(locus_link)
        # REmove all the files with that allele if present
        _ = [os.remove(x) for x in os.listdir('.') if locus['locus'] in x]
        
        with open(f"{locus['locus']}.fasta.gz", 'wb') as output_profile:
            output_profile.write(response_locus.read())
        lista_loci.append(f"{locus['locus']}")
        #slownik_alleli.update(get_feader(f"{locus['locus']}.fasta"))


    # donload the profile file
    _ = [os.remove(x) for x in os.listdir('.') if 'profiles.list' in x]
    profile_link = f"https://enterobase.warwick.ac.uk//schemes/{scheme_dir}/profiles.list.gz"
    response_profile = urlopen(profile_link)
    with open(f"profiles.list.gz", 'wb') as output_profile:
        output_profile.write(response_profile.read())
    execute_command(f"gunzip profiles.list.gz")
except HTTPError as Response_error:
    print(f"{Response_error.code} {Response_error.reason}. URL: {Response_error.geturl()}\\n Reason: {Response_error.read()}")

#np.save('slownik_wszystkich_alleli.npy', slownik_alleli, allow_pickle=True, fix_imports=True)

# indeks data with makeblastdb, more cpus the better
pool = Pool(cpus)

lista_indeksow = []
start = 0
step = len(lista_loci) // cpus

for i in range(cpus):
    end = start + step
    if i == (cpus - 1) or end > len(lista_loci):
        # if this is a last cpu or yoy passed the last element of a list, use as end len(lista_loci)
        end = len(lista_loci)
    lista_indeksow.append([start, end])
    start = end

jobs = []
for start, end in lista_indeksow:
    jobs.append(pool.apply_async(run_blast, (lista_loci, start, end)))
pool.close()
pool.join()

if not os.path.exists('local'):
    os.mkdir('local')

if not os.path.exists('local/profiles_local.list'):
    execute_command(f"head -1 profiles.list >> local/profiles_local.list")

