#!/usr/bin/env python3

import argparse
import csv
import json
import os
from warnings import warn


def parse_pangolin_output_csv2json(input: str, output: str, lan: str = 'en'):
    if not os.path.exists(input):
        msg = f"File {input} does not exist." if lan == 'en' else f"Plik {input} nie istnieje."
        output_json = {
            "status": "nie",
            "database_name": "Pangolin",
            "error_message": msg
        }
        with open(output, 'w') as f:
            f.write(json.dumps(output_json, ensure_ascii=False, indent=4))
    else:
        with open(input, 'r') as f:
            reader = csv.DictReader(f)
            data = list(reader)
        output_json = {
            "status": "tak",
            "database_name": "Pangolin",
            "database_version": data[0]['pangolin_version'],
            "program_name": "scorpion",
            "program_version": data[0]['scorpio_version'],
            "variant_name": data[0]['lineage'],
            "variant_qc_status": data[0]['qc_status'],
            "sequence_source": "full_genome"
        }
        with open(output, 'w') as f:
            f.write(json.dumps(output_json, ensure_ascii=False, indent=4))

def main():
    parser = argparse.ArgumentParser(description='Convert Pangolin output CSV to JSON')
    parser.add_argument('--input', type=str, required=True, help='Input CSV file')
    parser.add_argument('--output', type=str, required=True, help='Output JSON file')
    parser.add_argument('--lan', type=str, choices=['en', 'pl'], default='en',
                        help='Language for error messages')
    args = parser.parse_args()

    parse_pangolin_output_csv2json(args.input, args.output, args.lan)


if __name__ == '__main__':
    main()
