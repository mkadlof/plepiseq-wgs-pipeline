#!/usr/bin/env python3

"""
Download the latest Kraken2 database from Amazon S3.

Read more: https://benlangmead.github.io/aws-indexes/k2
"""

import argparse
import os
import os.path
import re
from datetime import datetime
from typing import List

import boto3
from botocore import UNSIGNED
from botocore.config import Config
import hashlib


# It employs the boto3 library to interact with the S3 service.
# It is assumed that the database name contains a date in the format YYYYMMDD.
# The script is using regex to extract the date from the database name, and then it finds the latest database.


def calculate_md5(file_path):
    """
    Calculate md5 sum of a file within python
    @param file_path:
    @type file_path:
    @return:
    @rtype:
    """
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def list_available_databases(bucket_name: str, prefix: str) -> List[str]:
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    try:
        response = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix)
        databases = [obj["Key"] for obj in response.get("Contents", [])]
        return databases
    except Exception as e:
        print(f"Error while listing databases: {e}")
        exit(1)


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
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    try:
        print(f"Downloading {file_name} from S3 to {local_path}")
        s3.download_file(bucket_name, file_name, local_path)

        real_path = os.path.dirname(local_path)
        os.system(f"tar -zxf {local_path} -C {real_path}")
        print(f"File downloaded successfully.")
        # Save md5 sum of tar.gz file to a file and remove the tar.gz file
        suma_md5 = calculate_md5(f"{local_path}")
        with open(f"{real_path}/current_md5.txt", "w") as f:
            f.write(f"{suma_md5}\n")
        os.system(f"rm {local_path}")

    except Exception as e:
        print(f"Error while downloading file: {e}")
        exit(1)

def check_updates(bucket_name: str, file_name: str, local_path: str):
    """
    Function checks if database version stored locally is identical to latest version available online
    @param bucket_name: Addres of database on AWS ?
    @type bucket_name: str
    @param file_name: This is actually a path to a .gz with database within "bucket"
    @type file_name: str
    @param local_path: This is a relative path where file from AWS is saved. Relative to the execution dir
    @type local_path: str
    @return: True if md5sums are different, Fasle  if latest version is available
    @rtype: bool
    """
    """Download a file from Amazon S3 and save it locally."""
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))

    # md5 files are stored in a differenct directory on aws so we need to adjust it
    # for
    # kraken/k2_standard_20240904.tar.gz
    # md5 is in
    # https://genome-idx.s3.amazonaws.com/kraken/standard_20240904/standard.md5
    real_name = os.path.basename(file_name)
    _ , typ, date = real_name.split('_')
    date = date.split('.')[0]

    # local_path = os.path.dirname(local_path)
    new_file_name = f"kraken/{typ}_{date}/{typ}.md5"
    print(f"Downloading md5 {new_file_name} to {local_path}/{typ}.md5")
    s3.download_file(bucket_name, new_file_name, f'{local_path}/{typ}.md5')

    with open(f"{local_path}/{typ}.md5") as f:
        for line in f:
            line = line.split()
            line[-1] = line[-1].rstrip()
            if "tar.gz" in line[-1]:
                novel_md5 = line[0]
                break

    old_md5=open(f"{local_path}/current_md5.txt").readlines()[0].rstrip()
    # Remove EVERYTHING from directory where kraken2 database will be saved
    os.system(f'rm {local_path}/{typ}.md5')
    # compare md5 sums of files return True if md5 sums are different and one need to perform update
    if novel_md5 == old_md5:
        return False
    else:
        return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("local_path", type=str, default=".")
    parser.add_argument("db_name", type=str,
                        choices=['standard', 'standard_08gb', 'standard_16gb', 'viral', 'minusb',
                                 'pluspf', 'pluspf_08gb', 'pluspf_16gb', 'pluspfp', 'pluspfp_08gb',
                                 'pluspfp_16gb', 'nt', 'eupathdb48'],
                        help="Database name (if unsure, use 'standard')")
    args = parser.parse_args()

    if os.path.dirname(args.local_path) == "" or os.path.dirname(args.local_path) == "/":
        print("Please provide full path as a first argument")
        exit(1)

    # Bucket name on Amazon S3
    bucket_name = "genome-idx"

    # Prefix for databases this is actually a directory where kraken2 database is stored on AWS
    prefix = "kraken/"

    # DB name regexp
    db_name_regexp = re.compile(prefix + r"k2_" + args.db_name + r"_(?P<date>\d{8})\.tar\.gz")

    # Get the list of available databases
    databases = list_available_databases(bucket_name, prefix)
    target_db = find_latest_database(databases, db_name_regexp)
    if not target_db:
        print("Database was not found.")
        exit(1)

    if not os.path.exists(args.local_path):
        print(f" Provided directory {args.local_path} does not exists. Sorry you must create it yourself")
        exit(1)


    # check if output directory contains two files current_md5.txt and database200mers.kmer_distrib
    # if yes, mostl likely we are trying to do an update

    if (not os.path.exists(f"{args.local_path}/current_md5.txt") and
            not os.path.exists(f"{args.local_path}/database200mers.kmer_distrib")):
        local_path = os.path.join(args.local_path, target_db.split("/")[-1])
        download_from_s3(bucket_name, target_db, local_path)
    else:
        print(f'Directory you provided {args.local_path} is not empty. Checking if new version of database is available')
        if check_updates(bucket_name, target_db, args.local_path):
            print('New version of database f{target_db} found, downloading new data')
            os.system(f'rm {args.local_path}/*')
            #  perform regular download
            local_path = os.path.join(args.local_path, target_db.split("/")[-1])
            download_from_s3(bucket_name, target_db, local_path)
        else:
            print(
                f"In direcotry: {args.local_path} latest version of database {target_db} is already present. Exiting")
            exit(1)


if __name__ == "__main__":
    main()
