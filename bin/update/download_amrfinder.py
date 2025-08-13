import os
import sys
import time
import logging
import subprocess
from datetime import datetime
from ftplib import FTP, error_perm
import io
import click

# --------------------
# Config & constants
# --------------------
MAX_RETRIES = 3
FTP_HOST = "ftp.ncbi.nlm.nih.gov"
FTP_DIR = "/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"
VERSION_FILE = "version.txt"
REQUIRED_FILES = [
    "AMR.LIB",
    "AMR_CDS.fa",
    "database_format_version.txt",
]

# --------------------
# Helpers
# --------------------
def execute_command(cmd: str) -> bool:
    # logging.info(f"Executing: {cmd}")
    try:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if stdout:
            logging.debug(stdout.decode(errors="replace"))
        if stderr:
            logging.debug(stderr.decode(errors="replace"))
        return proc.returncode == 0
    except Exception as e:
        logging.error(f"Command failed with exception: {e}")
        return False


def ftp_connect() -> FTP:
    ftp = FTP(FTP_HOST, timeout=60)
    ftp.login()
    ftp.cwd(FTP_DIR)
    return ftp


def ftp_read_text(ftp: FTP, filename: str) -> str:
    buf = io.BytesIO()
    ftp.retrbinary(f"RETR {filename}", buf.write)
    return buf.getvalue().decode("utf-8", errors="replace").strip()


def ftp_is_dir(ftp: FTP, name: str) -> bool:
    """Return True if 'name' on the current FTP path is a directory."""
    cur = ftp.pwd()
    try:
        ftp.cwd(name)
        ftp.cwd(cur)
        return True
    except error_perm:
        return False


def ftp_list_regular_files(ftp: FTP):
    try:
        names = ftp.nlst()
    except error_perm:
        names = []
    files = []
    for n in names:
        if n in (".", ".."):
            continue
        if ftp_is_dir(ftp, n):
            logging.info(f"Skipping directory on FTP: {n}")
            continue
        files.append(n)
    return files


def ftp_download_file(ftp: FTP, remote_name: str, local_path: str) -> None:
    with open(local_path, "wb") as fh:
        ftp.retrbinary(f"RETR {remote_name}", fh.write)


def is_fasta_and_dbtype(local_path: str):
    try:
        with open(local_path, "rt", errors="replace") as f:
            first = f.readline()
            if not first.startswith(">"):
                return (False, None)
            second = f.readline().strip()
            allowed = set("ATGCatgc")
            dbtype = "nucl" if set(second) <= allowed else "prot"
            return (True, dbtype)
    except Exception as e:
        logging.debug(f"Failed to inspect {local_path}: {e}")
        return (False, None)


def index_if_fasta(local_path: str) -> None:
    is_fa, dbtype = is_fasta_and_dbtype(local_path)
    if not is_fa:
        logging.debug(f"Not FASTA (skipping index): {os.path.basename(local_path)}")
        return
    execute_command(f"makeblastdb -in {local_path} -dbtype {dbtype}")


def required_files_present(output_dir: str) -> bool:
    return all(os.path.exists(os.path.join(output_dir, f)) for f in REQUIRED_FILES)


def clean_output_dir(output_dir: str) -> None:
    for name in os.listdir(output_dir):
        try:
            os.remove(os.path.join(output_dir, name))
        except Exception as e:
            logging.warning(f"Failed to remove {name}: {e}")


# --------------------
# CLI
# --------------------
@click.command()
@click.option('-o', '--output_dir', required=True, type=click.Path(),
              help='[REQUIRED] Directory where the AMRFinderPlus database will be stored')
@click.option('--force', is_flag=True, default=False,
              help='Remove all files in output_dir before downloading')
def main(output_dir: str, force: bool):
    os.makedirs(output_dir, exist_ok=True)

    log_file = os.path.join(output_dir, 'log.log')
    handlers = [
        logging.FileHandler(log_file, mode='w'),
        logging.StreamHandler(sys.stdout)
    ]
    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] %(levelname)s: %(message)s',
                        handlers=handlers)

    start = datetime.now()
    logging.info(f"Script started at {start.strftime('%Y-%m-%d %H:%M:%S')}")

    if force:
        logging.info("--force enabled: removing all files in output_dir")
        clean_output_dir(output_dir)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            logging.info(f"Attempt {attempt} of {MAX_RETRIES}...")
            ftp = ftp_connect()

            remote_version = ftp_read_text(ftp, VERSION_FILE)
            logging.info(f"Remote version: {remote_version}")

            local_version_path = os.path.join(output_dir, VERSION_FILE)
            local_version = None
            if os.path.exists(local_version_path):
                try:
                    with open(local_version_path, 'rt') as vf:
                        local_version = vf.read().strip()
                except Exception:
                    local_version = None

            # Decide whether to download
            if (not force) and local_version == remote_version and required_files_present(output_dir):
                logging.info(f"No update needed (local version {local_version} is up to date and required files present).")
                ftp.quit()
                break

            # Stamp remote version early (enables resume behavior)
            with open(local_version_path, 'wt') as vf:
                vf.write(remote_version + "")

            names = ftp_list_regular_files(ftp)
            if not names:
                raise RuntimeError("Could not list files on FTP server")

            for name in names:
                local_path = os.path.join(output_dir, os.path.basename(name))
                base = os.path.basename(name)
                # Remove previous copies (and index aux files) for this basename
                for fn in os.listdir(output_dir):
                    if fn.startswith(base):
                        try:
                            os.remove(os.path.join(output_dir, fn))
                        except Exception:
                            pass

                logging.info(f"Downloading: {name}")
                ftp_download_file(ftp, name, local_path)
                index_if_fasta(local_path)

            ftp.quit()

            # Validate presence of required files
            if not required_files_present(output_dir):
                raise RuntimeError("Required files missing after download")

            logging.info("Download and indexing completed successfully.")
            break

        except Exception as e:
            logging.warning(f"Attempt {attempt} failed: {e}")
            if attempt < MAX_RETRIES:
                logging.info("Cleaning output directory and retrying...")
                clean_output_dir(output_dir)
                time.sleep(120)
            else:
                logging.error("Failed to download AMRFinderPlus database after all retries.")
                sys.exit(1)

    end = datetime.now()
    logging.info(f"Script finished at {end.strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Total runtime: {str(end - start).split('.')[0]}")


if __name__ == '__main__':
    main()
