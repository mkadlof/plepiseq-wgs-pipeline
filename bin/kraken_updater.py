#!/usr/bin/env python3

"""
Download the latest Kraken2 database from Amazon S3.

Read more: https://benlangmead.github.io/aws-indexes/k2
"""

import argparse
import os.path
import re
from datetime import datetime
from typing import List

import boto3
from botocore import UNSIGNED
from botocore.config import Config


# It employs the boto3 library to interact with the S3 service.
# It is assumed that the database name contains a date in the format YYYYMMDD.
# The script is using regex to extract the date from the database name, and then it finds the latest database.

def list_available_databases(bucket_name: str, prefix: str) -> List[str]:
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    try:
        response = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix)
        databases = [obj["Key"] for obj in response.get("Contents", [])]
        return databases
    except Exception as e:
        print(f"Error while listing databases: {e}")
        return []


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
        print(f"File downloaded successfully.")
    except Exception as e:
        print(f"Error while downloading file: {e}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("db_name", type=str,
                        choices=['standard', 'standard_08gb', 'standard_16gb', 'viral', 'minusb',
                                 'pluspf', 'pluspf_08gb', 'pluspf_16gb', 'pluspfp', 'pluspfp_08gb',
                                 'pluspfp_16gb', 'nt', 'eupathdb48'],
                        help="Database name (if unsure, use 'standard')")
    parser.add_argument("local_path", type=str, default=".")
    args = parser.parse_args()

    # Bucket name on Amazon S3
    bucket_name = "genome-idx"

    # Prefix for databases
    prefix = "kraken/"

    # DB name regexp
    db_name_regexp = re.compile(prefix + r"k2_" + args.db_name + r"_(?P<date>\d{8})\.tar\.gz")

    # Get the list of available databases
    databases = list_available_databases(bucket_name, prefix)
    target_db = find_latest_database(databases, db_name_regexp)

    if target_db:
        local_path = os.path.join(args.local_path, target_db.split("/")[-1])
        # Check if the file already exists
        if not os.path.exists(local_path):
            download_from_s3(bucket_name, target_db, local_path)
        else:
            print(f"File {local_path} already exists. Skipping download.")
    else:
        print("Database was not found.")


if __name__ == "__main__":
    main()
