import requests
import subprocess
import sys
import time
import os

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


SPEC="pubmlst_campylobacter_seqdef"
DATABASE='C. jejuni / C. coli cgMLST v2'

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
    i = 0
    for line in profile.iter_lines():
        # Last column is header has an additional entry "LINcode" that should not be included
        if i == 0:
            line = list(map(lambda x: x.decode('utf-8', errors='replace'), line.split()))[:-1]
        else:
            line = list(map(lambda x: x.decode('utf-8', errors='replace'), line.split()))
        # replace "N" with 0, like in enterobase 
        line = ["0" if x == "N" else x for x in line]
        line = "\t".join(line)
        f.write(line + "\n")
        i = 1

#download and index fasta for each locus

loci = requests.get(scheme_link + '/loci')

for locus in loci.json()['loci']:
    locus_name = locus.split('/')[-1]
    print(f'Downloading data for locus: {locus_name}')
    fasta_file = requests.get(locus + '/alleles_fasta')
    with open(f'{locus_name}.fasta', 'w') as f:
        f.write(fasta_file.text)
    execute_command(f"makeblastdb -in {locus_name}.fasta -dbtype nucl")
    time.sleep(1)

# profilec ftom ch have "Ns" instead of negative values
execute_command("sed -i  s'/N/-1/'g profiles.list ")

if not os.path.exists('local'):
    os.mkdir('local')

if not os.path.exists('local/profiles_local.list'):
    execute_command(f"head -1 profiles.list >> local/profiles_local.list")

