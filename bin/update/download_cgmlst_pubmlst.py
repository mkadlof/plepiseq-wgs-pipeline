import requests
import subprocess
import sys
import time
import os
import logging
from multiprocessing import Pool
from datetime import datetime
import click

MAX_RETRIES = 3
REQUIRED_FILES = ["profiles.list", "CAMP1069.fasta", "CAMP0509.fasta"]

# ---- Helper Functions ----
def execute_command(command: str):
    logging.info(f"Executing command: {command}")
    try:
        result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = result.communicate()
        if stdout:
            logging.debug(stdout.decode())
        if stderr:
            logging.debug(stderr.decode())
        return result.returncode == 0
    except Exception as e:
        logging.error(f"Command execution failed: {e}")
        return False

def run_blast(lista_loci, start, end, output_dir):
    for locus in lista_loci[start:end]:
        execute_command(f"makeblastdb -in {os.path.join(output_dir, locus)}.fasta -dbtype nucl")
    return True

def is_data_complete(output_dir):
    for file in REQUIRED_FILES:
        if not os.path.exists(os.path.join(output_dir, file)):
            return False
    return True

# ---- Main CLI ----
@click.command()
@click.option('-c', '--cpus', default=4, show_default=True, help='Number of CPUs to use for BLAST indexing')
@click.option('-o', '--output_dir', help='[REQUIRED] Output directory', required=True, type=click.Path())
@click.option('-d', '--database', help='[REQUIRED] Genus-specific name of the database in Pubmlst ',
              type=click.Choice(["pubmlst_campylobacter_seqdef"]), required=True)
@click.option('-s', '--scheme_name', help='[REQUIRED] Name of the cgMLST scheme in Enterobase',
              type=click.Choice(["C. jejuni / C. coli cgMLST v2"]), required=True)
def main(cpus, output_dir, database, scheme_name):
    os.makedirs(output_dir, exist_ok=True)

    # Setup logging
    log_file_path = os.path.join(output_dir, "log.log")
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(log_file_path),
            logging.StreamHandler(sys.stdout)
        ]
    )

    start_time = datetime.now()
    logging.info(f"Script started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")


    if not execute_command('makeblastdb -version'):
        logging.error('makeblastdb was not found in PATH')
        sys.exit(1)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            logging.info(f"Download attempt {attempt}...")

            # find scheme link
            scheme_link = ''
            scheme_table = requests.get(f'https://rest.pubmlst.org/db/{database}/schemes')
            for scheme in scheme_table.json()['schemes']:
                if scheme_name == scheme['description']:
                    scheme_link = scheme['scheme']

            if not scheme_link:
                logging.error("Could not find scheme link")
                sys.exit(1)

            # download profiles
            logging.info("Downloading profiles")
            if os.path.exists(os.path.join(output_dir, 'profiles.list')):
                os.remove(os.path.join(output_dir, 'profiles.list'))

            profile = requests.get(scheme_link + '/profiles_csv')
            with open(os.path.join(output_dir, 'profiles.list'), 'w') as f:
                for line in profile.iter_lines():
                    line = list(map(lambda x: x.decode('utf-8', errors='replace'), line.split()))
                    line = ["0" if x == "N" else x for x in line]
                    f.write("\t".join(line[:1143]) + "\n")

            # download loci
            lista_loci = []
            loci = requests.get(scheme_link + '/loci')
            logging.info("Downloading loci")
            for i, locus in enumerate(loci.json()['loci']):

                if i % 100 == 0:
                    time.sleep(2)
                locus_name = locus.split('/')[-1]

                for f in os.listdir(output_dir):
                    if locus_name in f:
                        os.remove(os.path.join(output_dir, f))

                fasta_file = requests.get(locus + '/alleles_fasta')
                with open(os.path.join(output_dir, f'{locus_name}.fasta'), 'w') as f:
                    f.write(fasta_file.text)
                lista_loci.append(locus_name)

            if is_data_complete(output_dir):
                logging.info("All required files downloaded successfully.")
                break
            else:
                raise FileNotFoundError("Some required files are missing after download attempt.")

        except Exception as e:
            logging.warning(f"Attempt {attempt} failed: {str(e)}")
            if attempt < MAX_RETRIES:
                time.sleep(5)
            else:
                logging.error("Failed to download all required files after 3 attempts.")
                sys.exit(1)

    # Indexing
    logging.info("Starting indexing loci for BLAST...")
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

    local_dir = os.path.join(output_dir, 'local')
    os.makedirs(local_dir, exist_ok=True)

    if not os.path.exists(os.path.join(local_dir, 'profiles_local.list')):
        execute_command(f"head -1 {os.path.join(output_dir, 'profiles.list')} >> {os.path.join(local_dir, 'profiles_local.list')}")

    end_time = datetime.now()
    logging.info(f"Script finished at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Total runtime: {str(end_time - start_time).split('.')[0]}")

if __name__ == '__main__':
    main()
