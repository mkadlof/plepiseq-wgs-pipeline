import requests
import subprocess
import sys
import time
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

def run_blast(lista_loci, start, end):
    # index the file
    for locus in lista_loci[start:end]:
        execute_command(f"makeblastdb -in {locus}.fasta -dbtype nucl")
    return True


SPEC="pubmlst_campylobacter_seqdef"
DATABASE='C. jejuni / C. coli cgMLST v2'
cpus = int(sys.argv[1])

# check if makeblastd is in path

if not execute_command('makeblastdb  -version'):
    print('makeblastdb was not foundh')
    sys.exit(1)

# find scheme link
scheme_link = ''
scheme_table = requests.get('https://rest.pubmlst.org/db/' f'{SPEC}' + '/schemes')
for scheme in scheme_table.json()['schemes']:
    if DATABASE == scheme['description']:
        scheme_link =  scheme['scheme']
        
# extract profile
print(f'Downloading profiles')
profile = requests.get(scheme_link + '/profiles_csv')
with open('profiles.list', 'w') as f:
    for line in profile.iter_lines():
        # Last column is LINcode that should not be included
        line = list(map(lambda x: x.decode('utf-8', errors='replace'), line.split()))
        # replace "Ns" with 0
        line = ["0" if x == "N" else x for x in line]
        line = "\t".join(line[:1143])
        f.write(line + "\n")

#download and index fasta for each locus

loci = requests.get(scheme_link + '/loci')

lista_loci = []

i = 0
for locus in loci.json()['loci']:
    if i % 100 == 0:
        time.sleep(1)
    locus_name = locus.split('/')[-1]
    print(f'Downloading data for locus: {locus_name}')
    lista_loci.append(locus_name)
    fasta_file = requests.get(locus + '/alleles_fasta')
    with open(f'{locus_name}.fasta', 'w') as f:
        f.write(fasta_file.text)
    i += 1

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

