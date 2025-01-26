from urllib.request import urlopen
from urllib.error import HTTPError
import urllib
import base64
import json
import os
import subprocess
import time
from Bio import SeqIO
import numpy as np


API_TOKEN=open('/home/update/enterobase_api.txt').readlines()[0].rstrip()
DATABASE="senterica" # ecoli, yersinia, mcatarrhalis API 2.0
scheme_name="MLST_Achtman" # according to API 2.0 
scheme_dir="Salmonella.Achtman7GeneMLST" # take a look inside https://enterobase.warwick.ac.uk//schemes/ for names of specific schemes


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

if os.path.exists('slownik_wszystkich_alleli.npy'):
    # slownik_wszystkich_alleli.npy trzyma informacje jakie allele i ich wersje sa w bazie struktura to {"NAZWA_ALLELU':[1,2,3,4,5]; "NAZWA_KOLEJNEGO_ALLELU: [1,2] }
    # idea jest taka jesli jest ten slownik tzn ze pobrano juz chociaz raz allele wiec sa tez pliki fasta
    # potem poprzez interfase w /api/v2.0/{database}/{scheme}/alleles 
    # pobieramy liste wszystkich alleli , sprawdzamy jakei sa wersje allelu i nowe sekwencje appendujemy do plikow fasta ( a potem indeksujemy makeblastdb, po usunieciu starcych indeksow)
    pass
    # Update  jest jednak bez sensu bo sprowadza sie do odpytania bazy sprawdzenia co jest nowe i sciagniecia sekwencji
    # Kiedy mozna po prostu sciagnac sekwencje 

    # Nie ma slownika_alleli pobieramy wszystkie allele "na pale" z odpowiedniej bazy z https://enterobase.warwick.ac.uk/schemes/
    # tworzymy tez slownik i zapisujemy jako slownik_wszystkich_alleli.npy

address = f'https://enterobase.warwick.ac.uk/api/v2.0/{DATABASE}/{scheme_name}/loci?limit=10000&scheme={scheme_name}&offset=0'
slownik_alleli = {}
try:
    os.remove('all_allels.fasta')
    os.remove("MLST_Achtman_ref.fasta")
except FileNotFoundError:
    pass

try:
    response = urlopen(__create_request(address))
    data = json.load(response)
    for locus in data['loci']:
        time.sleep(2)
        print(locus['locus'])
        locus_link = f"https://enterobase.warwick.ac.uk//schemes/{scheme_dir}/{locus['locus']}.fasta.gz"
            # download the locus_address
        response_locus = urlopen(locus_link)
        # REmove all the files with that allele if present
        _ = [os.remove(x) for x in os.listdir('.') if locus['locus'] in x]
        
        with open(f"{locus['locus']}.fasta.gz", 'wb') as output_profile:
            output_profile.write(response_locus.read())
        # index the file
        execute_command(f"gunzip {locus['locus']}.fasta.gz")
        execute_command(f"makeblastdb -in {locus['locus']}.fasta -dbtype nucl")
        slownik_alleli.update(get_feader(f"{locus['locus']}.fasta"))
    # donload the profile file
    _ = [os.remove(x) for x in os.listdir('.') if 'profiles.list' in x]
    profile_link = f"https://enterobase.warwick.ac.uk//schemes/{scheme_dir}/profiles.list.gz"
    response_profile = urlopen(profile_link)
    with open(f"profiles.list.gz", 'wb') as output_profile:
        output_profile.write(response_profile.read())
    execute_command(f"gunzip profiles.list.gz")

    # merge all fasta into one file
    execute_command(f"cat *fasta > all_allels.fasta")
    
    # download reference fasta
    # File is created when preparing MLSTdb 
    #fasta_link = f"https://enterobase.warwick.ac.uk//schemes/{scheme_dir}/MLST_Achtman_ref.fasta"
    #fasta_response = urlopen(fasta_link)
    #with open(f"MLST_Achtman_ref.fasta", 'wb') as output_profile:
    #    output_profile.write(fasta_response.read())

except HTTPError as Response_error:
    print(f"{Response_error.code} {Response_error.reason}. URL: {Response_error.geturl()}\\n Reason: {Response_error.read()}")


np.save('slownik_wszystkich_alleli.npy', slownik_alleli, allow_pickle=True, fix_imports=True)
if not os.path.exists('local'):
    os.mkdir('local')

if not os.path.exists('local/profiles_local.list'):
    execute_command(f"head -1 profiles.list >> local/profiles_local.list")
