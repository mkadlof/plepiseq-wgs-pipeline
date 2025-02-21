#!/usr/bin/env python3
"""
Skrypt do zwracania statusu w module run_initial_mlst*
"""

import sys
import click
import json


@click.command()
@click.option('-i', '--input_file', help='[INPUT] a path to a file MLST prediction ',
              type=click.Path(), default="fastqc_data.txt", required = False )
@click.option('-x', '--min_number', help='[INPUT] QC parameter minimal number of unique loci',
              type=int,  required = False )
@click.option('-s', '--status', help='[INPUT] PREDEFINED status that is transferred to an output json. '
                                     'If this status was either nie or blad fastqc will not run',
              type=click.Choice(['tak', 'nie', 'blad'], case_sensitive=False),  required=True)
@click.option('-r', '--error', help='[INPUT] PREDEFINED error message that is put in json. '
                                    'Only used when status was set to nie or blad',
              type=str,  required=False, default="")
@click.option('-o', '--output', help='[Output] Name of a file with json output',
              type=str,  required=False)
def main_program(input_file, min_number, status, output="", error=""):
    slownik_loci = {}  # slownik zawiera jako klucz naze loci, a jako wartosc ilosc alleli
    # ze 100% seq identity i 100% pokryciem
    if status != "tak":
        json_dict = {"scheme_name": "MLST_cge",
                     "status": status,
                     "error_message": error}
        status = "nie" # status przekazany do modulow ponizej
    else:
        with open(input_file) as f:
            for line in f:
                line = line.split("\t")
                if "#" in line[0]:
                    continue
                loci_nazwa = line[0].split("_")[0]
                allel_seq_id = float(line[4])
                allel_coverage = float(line[5])
                loci_coverage = float(line[6])
                loci_seq_id = float(line[7])
                if loci_nazwa not in slownik_loci.keys():
                    slownik_loci[loci_nazwa] = 0

                if allel_seq_id == 100 and allel_coverage == 100 and loci_coverage == 100 and loci_seq_id == 100:
                    slownik_loci[loci_nazwa] += 1

        unikalne_loci = len([True for x, y in slownik_loci.items() if y == 1])
        duplicated_loci = len([True for x, y in slownik_loci.items() if y > 1])
        missing_loci = 7 - unikalne_loci - duplicated_loci

        if unikalne_loci >= float(min_number):
            status = "tak"
            json_dict = {"scheme_name": "MLST_cge",
                         "status": status,
                         "missing_allels_value": missing_loci,
                         "duplicated_allels_value": duplicated_loci,
                         "unique_allels_value": unikalne_loci}

        else:
            status = "blad" # status do jsona
            error = f"The number of unique loci in sample is {unikalne_loci}, which is below designed threshold"
            json_dict = {"scheme_name": "MLST_cge",
                         "status": status,
                         "error_message": error}
            status = "nie"# status przekazywany modulow nizej

    with open(output, 'w') as f1:
        f1.write(json.dumps(json_dict, indent = 4))
    print(status)
    return status


if __name__ == '__main__':
    # The main program returns 3 variables: "status", number_of_reads in a sample, and median_quality_of_reads
    # These can be used to determine QC status in tha module
    if len(sys.argv) == 1:
        main_program(['--help'])
    else:
        main_program(sys.argv[1:])
