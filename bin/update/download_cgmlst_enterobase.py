import urllib.request
from urllib.error import HTTPError
import base64
import json
import os
import subprocess
import time
from Bio import SeqIO
import sys
from multiprocessing import Pool
import logging
from datetime import datetime

# ---- Logging setup ----
logging.basicConfig(
    format='[%(asctime)s] %(levelname)s: %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S'
)

def __create_request(request_str):
    base64string = base64.b64encode(f'{API_TOKEN}: '.encode('utf-8'))
    headers = {"Authorization": f"Basic {base64string.decode()}"}
    request = urllib.request.Request(request_str, None, headers)
    return request

def execute_command(command: str):
    """Run a shell command and wait for it to finish."""
    logging.info(f"Executing command: {command}")
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stdout:
        logging.debug(stdout.decode())
    if stderr:
        logging.debug(stderr.decode())
    return process.returncode == 0

def run_blast(lista_loci, start, end):
    for locus in lista_loci[start:end]:
        execute_command(f"gunzip {locus}.fasta.gz")
        execute_command(f"makeblastdb -in {locus}.fasta -dbtype nucl")
    return True

# ---- Start of script ----
start_time = datetime.now()
logging.info(f"Script started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

API_TOKEN = open('/home/update/enterobase_api.txt').readlines()[0].rstrip()
DATABASE = sys.argv[1]        # ecoli, yersinia, etc.
scheme_name = sys.argv[2]     # e.g., Achtman7
scheme_dir = sys.argv[3]      # e.g., Achtman7
cpus = int(sys.argv[4])       # number of CPUs

address = f'https://enterobase.warwick.ac.uk/api/v2.0/{DATABASE}/{scheme_name}/loci?limit=10000&scheme={scheme_name}&offset=0'
lista_loci = []

try:
    logging.info(f"Requesting loci list from EnteroBase database: {DATABASE} ...")
    response = urllib.request.urlopen(__create_request(address))
    data = json.load(response)

    i = 0
    for locus in data['loci']:
        if i % 100 == 0:
            time.sleep(1)
        i += 1

        locus_link = f"https://enterobase.warwick.ac.uk//schemes/{scheme_dir}/{locus['locus']}.fasta.gz"
        response_locus = urllib.request.urlopen(locus_link)

        # Remove existing files matching the locus
        _ = [os.remove(x) for x in os.listdir('.') if locus['locus'] in x]

        with open(f"{locus['locus']}.fasta.gz", 'wb') as f_out:
            f_out.write(response_locus.read())
        lista_loci.append(f"{locus['locus']}")

    # Download profile file
    logging.info("Downloading profile file...")
    _ = [os.remove(x) for x in os.listdir('.') if 'profiles.list' in x]
    profile_link = f"https://enterobase.warwick.ac.uk//schemes/{scheme_dir}/profiles.list.gz"
    response_profile = urllib.request.urlopen(profile_link)
    with open("profiles.list.gz", 'wb') as f_out:
        f_out.write(response_profile.read())
    execute_command("gunzip profiles.list.gz")

except HTTPError as e:
    logging.error(f"{e.code} {e.reason}. URL: {e.geturl()} Reason: {e.read()}")
    sys.exit(1)

# ---- Indexing for BLAST ----
logging.info("Starting indexing loci for BLAST...")
pool = Pool(cpus)
lista_indeksow = []
start = 0
step = len(lista_loci) // cpus

for i in range(cpus):
    end = start + step
    if i == (cpus - 1) or end > len(lista_loci):
        end = len(lista_loci)
    lista_indeksow.append([start, end])
    start = end

jobs = []
for start, end in lista_indeksow:
    jobs.append(pool.apply_async(run_blast, (lista_loci, start, end)))

pool.close()
pool.join()

# ---- Final steps ----
if not os.path.exists('local'):
    os.mkdir('local')

if not os.path.exists('local/profiles_local.list'):
    execute_command("head -1 profiles.list >> local/profiles_local.list")

end_time = datetime.now()
logging.info(f"Script finished at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
logging.info(f"Total runtime: {str(end_time - start_time).split('.')[0]}")

