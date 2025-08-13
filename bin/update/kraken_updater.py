#!/usr/bin/env python3

"""
Download the latest Kraken2 database from Amazon S3.

Read more: https://benlangmead.github.io/aws-indexes/k2
"""

import os
import os.path
import re
import sys
from datetime import datetime
from typing import List
import hashlib
import logging
import time

import boto3
from botocore import UNSIGNED
from botocore.config import Config
from botocore.exceptions import BotoCoreError, ClientError

import click

# --------------------
# Constants
# --------------------
DB_NAMES = [
    'standard', 'standard_08gb', 'standard_16gb', 'viral', 'minusb',
    'pluspf', 'pluspf_08gb', 'pluspf_16gb', 'pluspfp', 'pluspfp_08gb',
    'pluspfp_16gb', 'nt', 'eupathdb48'
]

S3_BUCKET = "genome-idx"
S3_PREFIX = "kraken/"

# --------------------
# Logging
# --------------------

def setup_logging(output_dir: str) -> None:
    """
    Simple function used to config logging
    """
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, 'log.log')
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(log_path, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )

# --------------------
# Helpers
# --------------------

def s3_client_unsigned():
    """
    We log to boto3 multiple times
    """
    return boto3.client("s3", config=Config(signature_version=UNSIGNED))


def check_s3_connectivity(bucket: str, prefix: str, attempts: int = 3, interval_sec: int = 120) -> None:
    """Attempt to list one object from S3 up to `attempts` times.
    Exits the program with a non-zero status if all attempts fail.
    """
    cli = s3_client_unsigned()
    for attempt in range(1, attempts + 1):
        try:
            logging.info(f"Checking S3 connectivity (attempt {attempt}/{attempts})...")
            resp = cli.list_objects_v2(Bucket=bucket, Prefix=prefix, MaxKeys=1)
            logging.info("S3 connectivity confirmed.")
            return
        except (BotoCoreError, ClientError) as e:
            logging.warning(f"S3 connectivity check failed: {e}")
            if attempt < attempts:
                logging.info(f"Retrying in {interval_sec} seconds...")
                time.sleep(interval_sec)
            else:
                logging.error("Unable to establish S3 connectivity after multiple attempts. Exiting.")
                sys.exit(1)


def calculate_md5(file_path):
    """
    Calculate md5 sum of a file within python
    """
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def list_available_databases(bucket_name: str, prefix: str) -> List[str]:
    cli = s3_client_unsigned()
    try:
        paginator = cli.get_paginator("list_objects_v2")
        keys: list[str] = []
        for page in paginator.paginate(Bucket=bucket_name, Prefix=prefix):
            for obj in page.get("Contents", []):
                keys.append(obj["Key"])
        return keys
    except Exception as e:
        logging.error(f"Error while listing databases: {e}")
        sys.exit(1)


def find_latest_database(databases: List[str], db_name_regexp: re.Pattern) -> str:
    """Find the latest database in the list of databases.
    It is assumed that the database name contains a date in the format YYYYMMDD."""
    latest_date = None
    latest_database = None
    for database in databases:
        match = db_name_regexp.match(database)
        if match:
            date_str = match.group("date")
            date = datetime.strptime(date_str, "%Y%m%d")
            if latest_date is None or date > latest_date:
                latest_date = date
                latest_database = database
    return latest_database


def download_from_s3(bucket_name: str, file_name: str, local_path: str):
    """Download a file from Amazon S3 and save it locally."""
    s3 = s3_client_unsigned()
    try:
        logging.info(f"Downloading {file_name} from S3 to {local_path}")
        s3.download_file(bucket_name, file_name, local_path)

        real_path = os.path.dirname(local_path)
        os.system(f"tar -zxf {local_path} -C {real_path}")
        logging.info("File downloaded and extracted successfully.")
        # Save md5 sum of tar.gz file to a file and remove the tar.gz file
        suma_md5 = calculate_md5(f"{local_path}")
        with open(f"{real_path}/current_md5.txt", "w") as f:
            f.write(f"{suma_md5}\n")
        os.system(f"rm {local_path}")

    except Exception as e:
        logging.error(f"Error while downloading file: {e}")
        sys.exit(1)


def check_updates(bucket_name: str, file_name: str, local_path: str):
    """
    Function checks if database version stored locally is identical to latest version available online.

    @return: True if md5sums are different (update needed), False if latest version is already present.
    """
    s3 = s3_client_unsigned()

    # Example:
    # kraken/k2_standard_20240904.tar.gz  -> md5 at  kraken/standard_20240904/standard.md5
    real_name = os.path.basename(file_name)
    _, typ, date = real_name.split('_')
    date = date.split('.')[0]

    new_file_name = f"kraken/{typ}_{date}/{typ}.md5"
    md5_tmp = os.path.join(local_path, f"{typ}.md5")
    logging.info(f"Downloading md5 {new_file_name} to {md5_tmp}")
    s3.download_file(bucket_name, new_file_name, md5_tmp)

    novel_md5 = None
    with open(md5_tmp) as f:
        for line in f:
            line = line.split()
            line[-1] = line[-1].rstrip()
            if "tar.gz" in line[-1]:
                novel_md5 = line[0]
                break

    try:
        old_md5 = open(f"{local_path}/current_md5.txt").readlines()[0].rstrip()
    except FileNotFoundError:
        # No current_md5.txt file, force update
        logging.info("No current_md5.txt found; will download latest database.")
        os.remove(md5_tmp)
        return True

    # Remove the temporary md5 file
    os.remove(md5_tmp)

    # compare md5 sums of files; return True if they differ and we need to update
    return (novel_md5 or "") != old_md5


@click.command()
@click.option(
    "-o", "--local_path",
    required=True,
    type=click.Path(),
    help="[REQUIRED] Full path to the destination directory",
)
@click.option(
    "-d", "--db_name",
    required=True,
    type=click.Choice(DB_NAMES),
    help="Database name (if unsure, use 'standard')",
)
def main(local_path: str, db_name: str):
    # Logging setup
    setup_logging(local_path)
    logging.info("Starting Kraken2 DB updater")

    # Argument validation (kept same spirit as original)
    if os.path.dirname(local_path) == "" or os.path.dirname(local_path) == "/":
        logging.error("Please provide full path as a first argument")
        sys.exit(1)

    # Check connection to S3
    check_s3_connectivity(S3_BUCKET, S3_PREFIX, attempts=3, interval_sec=120)

    # Bucket name on Amazon S3
    bucket_name = S3_BUCKET

    # Prefix for databases (S3 directory where Kraken2 db is stored)
    prefix = S3_PREFIX

    # DB name regexp
    db_name_regexp = re.compile(prefix + r"k2_" + db_name + r"_(?P<date>\d{8})\.tar\.gz")

    # Get the list of available databases
    databases = list_available_databases(bucket_name, prefix)

    target_db = find_latest_database(databases, db_name_regexp)

    if not target_db:
        logging.error("Database was not found.")
        sys.exit(1)

    if not os.path.exists(local_path):
        logging.error(f"Provided directory {local_path} does not exist. Please create it first.")
        sys.exit(1)

    # If empty (no current_md5.txt and no database200mers.kmer_distrib), treat as fresh download
    if (not os.path.exists(f"{local_path}/current_md5.txt") and
            not os.path.exists(f"{local_path}/database200mers.kmer_distrib")):
        dest_tar = os.path.join(local_path, target_db.split("/")[-1])
        download_from_s3(bucket_name, target_db, dest_tar)
    else:
        logging.info(f"Directory you provided {local_path} is not empty. Checking if new version of database is available")
        if check_updates(bucket_name, target_db, local_path):
            logging.info(f"New version of database {target_db} found, downloading new data")
            os.system(f'rm {local_path}/*')
            dest_tar = os.path.join(local_path, target_db.split("/")[-1])
            download_from_s3(bucket_name, target_db, dest_tar)
        else:
            logging.info(f"In directory: {local_path} latest version of database {target_db} is already present. Exiting")
            sys.exit(0)


if __name__ == "__main__":
    main()
