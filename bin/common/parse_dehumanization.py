#!/usr/bin/env python3

"""
This script parse the results of picard_wgsmetrics module in the viral pipeline to
produce a valid json (the viral_genom_data group of fileds)
"""

import sys
import click
import json


@click.command()
@click.option('-i', '--input_fastas_list', help='[INPUT] a file with a list of fastas produced by this module ',
              type=click.Path(), required=False)
@click.option('-k', '--output_path', help='[INPUT] a path where all the relevant files are put by nextflow',
              type=str)
@click.option('-s', '--status', help='[INPUT] PREDEFINED status that is transferred to an output json. '
                                     'If this status was either nie or blad fastqc will not run',
              type=click.Choice(['tak', 'nie', 'blad'], case_sensitive=False),  required=True)
@click.option('-r', '--error', help='[INPUT] PREDEFINED error message that is put in json. '
                                    'Only used when status was set to nie or blad',
              type=str,  required=False, default="")
@click.option('-o', '--output', help='[Output] Name of a file with json output',
              type=str,  required=True)
def main_program(status, output, input_fastas_list, output_path, error=""):
    total_length_value = 0
    number_of_Ns_value = 0
    if status != "tak":
        json_output = {"status": status,
                       "error_message": error}
    else:
        file_data = []
        with open(input_fastas_list) as f:
            file_data = [f"{output_path}/{x.rstrip()}" for x in f]

        json_output = {
            "status": status,
            "dehumanized_fastq_file_list" : file_data
        }
    with open(output, 'w') as f1:
        f1.write(json.dumps(json_output, indent=4))

    return True


if __name__ == '__main__':
    if len(sys.argv) == 1:
        main_program(['--help'])
    else:
        main_program(sys.argv[1:])
