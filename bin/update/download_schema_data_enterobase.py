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
import click
from Bio import SeqIO


# Here we define numbr of retries and expected files
MAX_RETRIES = 3
REQUIRED_FILES_CGMLST = {
    "senterica": ["STMMW_17971.fasta.gz", "t1733.fasta.gz"],
    "ecoli": ["b0784.fasta.gz", "NCTC12130_00627.fasta.gz"]
}

REQUIRED_FILES_MLST = {
    "senterica": ["aroC.fasta.gz", "dnaN.fasta.gz", "purE.fasta.gz"],
    "ecoli": ["adk.fasta.gz", "fumC.fasta.gz", "recA.fasta.gz"]
}

# ---- Helper Functions ----
def __create_request(request_str, api_token):
    base64string = base64.b64encode(f'{api_token}: '.encode('utf-8'))
    headers = {"Authorization": f"Basic {base64string.decode()}"}
    return urllib.request.Request(request_str, None, headers)

def execute_command(command: str):
    """
    generic function to execute a bash command
    """
    # logging.info(f"Executing command: {command}")
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stdout:
        logging.debug(stdout.decode())
    if stderr:
        logging.debug(stderr.decode())
    return process.returncode == 0

def run_blast(lista_loci, start, end, output_dir):
    """
    run blast
    """
    for locus in lista_loci[start:end]:
        execute_command(f"gunzip {os.path.join(output_dir, locus)}.fasta.gz")
        execute_command(f"makeblastdb -in {os.path.join(output_dir, locus)}.fasta -dbtype nucl")
    return True

def is_data_complete(database: str, output_dir: str, scheme_name:str):
    """
    Check if output directory contains relevant files
    """
    if not os.path.exists(os.path.join(output_dir, "profiles.list")):
        return False

    # Check expected files given organism and schema type
    if scheme_name in ["cgMLST_v2", "cgMLST"]:
        required = REQUIRED_FILES_CGMLST.get(database, [])
    elif scheme_name == "MLST_Achtman":
        required = REQUIRED_FILES_MLST.get(database, [])

    for file in required:
        if not os.path.exists(os.path.join(output_dir, file)):
            return False
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

# ---- Main Execution ----
@click.command()
@click.option('-d', '--database', help='[REQUIRED] Genus-specific name of the database in Enterobase ',
              type=click.Choice(['senterica', 'ecoli']), required=True)
@click.option('-s', '--scheme_name', help='[REQUIRED] Name of the cgMLST scheme in Enterobase',
              type=click.Choice(['cgMLST_v2', 'cgMLST', 'MLST_Achtman']), required=True)
@click.option('-r', '--scheme_dir', help='[REQUIRED] Schema directory in EnteroBase',
              type=click.Choice(['Salmonella.cgMLSTv2', 'Escherichia.cgMLSTv1', 'Escherichia.Achtman7GeneMLST', 'Salmonella.Achtman7GeneMLST']), required=True)
@click.option('-c', '--cpus', default=4, show_default=True, help='Number of CPUs to use for BLAST indexing')
@click.option('-t', '--api_token_file',
              default='/home/update/enterobase_api.txt', show_default=True,
              type=click.Path(), help='Path to EnteroBase API token file')
@click.option('-o', '--output_dir', help='[REQUIRED] output directory', required=True, type=click.Path())
def main(database, scheme_name, scheme_dir, cpus, api_token_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Setup logging to file in output directory
    log_file_path = os.path.join(output_dir, "log.log")
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(log_file_path, mode = 'w'),
            logging.StreamHandler(sys.stdout)
        ]
    )

    if scheme_name == "MLST_Achtman":
        execute_command('rm all_allels.fasta || true')
        execute_command('rm MLST_Achtman_ref.fasta || true')


    start_time = datetime.now()
    logging.info(f"Script started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    with open(api_token_file) as f:
        api_token = f.readline().strip()

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            logging.info(f"Download attempt {attempt}...")
            lista_loci = []

            address = f'https://enterobase.warwick.ac.uk/api/v2.0/{database}/{scheme_name}/loci?limit=10000&scheme={scheme_name}&offset=0'

            response = urllib.request.urlopen(__create_request(address, api_token))
            data = json.load(response)

            logging.info("Downloading loci files...")
            for i, locus in enumerate(data['loci']):
                if i % 100 == 0:
                    # In some cases API doesn't like to be queried to many times
                    time.sleep(2)
                locus_link = f"https://enterobase.warwick.ac.uk//schemes/{scheme_dir}/{locus['locus']}.fasta.gz"
                response_locus = urllib.request.urlopen(locus_link)

                for f in os.listdir(output_dir):
                    if locus['locus'] in f:
                        os.remove(os.path.join(output_dir, f))
                with open(os.path.join(output_dir, f"{locus['locus']}.fasta.gz"), 'wb') as f_out:
                    f_out.write(response_locus.read())
                lista_loci.append(f"{locus['locus']}")

            logging.info("Downloading profile file...")
            for f in os.listdir(output_dir):
                if 'profiles.list' in f:
                    os.remove(os.path.join(output_dir, f))

            profile_link = f"https://enterobase.warwick.ac.uk//schemes/{scheme_dir}/profiles.list.gz"
            response_profile = urllib.request.urlopen(profile_link)
            with open(os.path.join(output_dir, "profiles.list.gz"), 'wb') as f_out:
                f_out.write(response_profile.read())
            execute_command(f"gunzip {os.path.join(output_dir, 'profiles.list.gz')}")

            if is_data_complete(database, output_dir, scheme_name):
                logging.info("All required files downloaded successfully.")
                break
            else:
                raise FileNotFoundError("Some required files are missing after download attempt.")

        except Exception as e:
            logging.warning(f"Attempt {attempt} failed: {str(e)}")
            if attempt < MAX_RETRIES:
                time.sleep(30)
            else:
                logging.error("Failed to download all required files after 3 attempts.")
                sys.exit(1)

    logging.info("Starting indexing loci for BLAST...")

    # Turn off multiprocessing for 7-gene MLST
    if scheme_name == "MLST_Achtman":
        cpus = 1

    pool = Pool(cpus)
    lista_indeksow = []
    step = len(lista_loci) // cpus
    start = 0

    for i in range(cpus):
        end = start + step
        if i == (cpus - 1) or end > len(lista_loci):
            end = len(lista_loci)
        lista_indeksow.append((start, end))
        start = end

    jobs = [pool.apply_async(run_blast, (lista_loci, s, e, output_dir)) for s, e in lista_indeksow]
    pool.close()
    pool.join()

    if scheme_name == "MLST_Achtman":
        execute_command(f"cat *fasta > all_allels.fasta")

    # ---- Final steps ----
    # In case we run our script on an empty directory we must set up files with a local 'database'
    local_dir = os.path.join(output_dir, 'local')
    os.makedirs(local_dir, exist_ok=True)

    if not os.path.exists(os.path.join(local_dir, 'profiles_local.list')):
        execute_command(f"head -1 {os.path.join(output_dir, 'profiles.list')} >> {os.path.join(local_dir, 'profiles_local.list')}")

    end_time = datetime.now()
    logging.info(f"Script finished at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Total runtime: {str(end_time - start_time).split('.')[0]}" )

if __name__ == '__main__':
    main()
