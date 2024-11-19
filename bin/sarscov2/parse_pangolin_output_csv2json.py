#!/usr/bin/env python3

import argparse
import csv
import json
import os
from warnings import warn


def parse_pangolin_output_csv2json(input: str, output: str):
    if not os.path.exists(input):
        warn(f"File {input} does not exist.")
        output_json = {"status": "nie",
                       "database_name": "Pangolin",
                       "error_message": f"File {input} does not exist."}
        with open(output, 'w') as f:
            json.dump(output_json, f)
    else:
        with open(input, 'r') as f:
            reader = csv.DictReader(f)
            data = list(reader)
        output_json = {"status": "tak",
                       "database_name": "Pangolin",
                       "database_version": data[0]['pangolin_version'],
                       "program_name": "scorpion",
                       "program_version": data[0]['scorpio_version'],
                       "variant_name": data[0]['lineage'],
                       "variant_qc_status": data[0]['qc_status'],
                       "sequence_source": "full_genome",
                       }
        with open(output, 'w') as f:
            json.dump(output_json, f, indent=4)

    print(f"File {output} saved.")


def main():
    parser = argparse.ArgumentParser(description='Convert Pangolin output CSV to JSON')
    parser.add_argument('input', type=str, help='Input CSV file')
    parser.add_argument('output', type=str, help='Output JSON file')
    args = parser.parse_args()

    parse_pangolin_output_csv2json(args.input, args.output)


if __name__ == '__main__':
    main()
