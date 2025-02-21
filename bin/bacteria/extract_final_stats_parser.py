#!/usr/bin/env python3
"""
Skrypt do zwracania statusu i json w module extract_final_stats
"""

import sys
import click
import json


@click.command()
@click.option('-i', '--input_file_filtered', help='[INPUT] a path to an input file with summary '
                                                  'statistics of a genome that is based on filtered contigs ',
              type=click.Path(), default="fastqc_data.txt", required=False)
@click.option('-j', '--input_file_unfiltered', help='[INPUT] a path to an input file with summary '
                                                    'statistics of a genome that is based on unfiltered contigs ',
              type=click.Path(), default="fastqc_data.txt", required=False)
@click.option('-l', '--nvalue', help='[INPUT] QC parameter minimal value of N50',
              type=int,  required=True)
@click.option('-n', '--max_contigs', help='[INPUT] QC parameter maximal number of contigs',
              type=int,  required=True)
@click.option('-g', '--genome_length', help='[INPUT] Length of an expected genome',
              type=int,  required=True)
@click.option('-c', '--completness', help='[INPUT] QC parameter. Fraction of genome completness',
              type=float,  required=True)
@click.option('-p', '--coverage', help='[INPUT] QC parameter. Minimal average coverage',
              type=int,  required=True)
@click.option('-s', '--status', help='[INPUT] PREDEFINED status that is transferred to an output json. '
                                     'If this status was either nie or blad fastqc will not run',
              type=click.Choice(['tak', 'nie', 'blad'], case_sensitive=False),  required=True)
@click.option('-r', '--error', help='[INPUT] PREDEFINED error message that is put in json. '
                                    'Only used when status was set to nie or blad',
              type=str,  required=False, default="")
@click.option('-o', '--output', help='[Output] Name of a file with json output',
              type=str,  required=False)
def main_program(input_file_filtered, input_file_unfiltered, max_contigs, nvalue, genome_length,
                 completness, coverage, status, output="", error=""):
    if status != "tak":
        json_output = [{"step_name": "pre-filtering",
                        "status": status,
                        "error_message": error},
                       {"step_name": "post-filtering",
                        "status": status,
                        "error_message": error}]
    else:
        json_output = []

        # initial values of all parameters passed to a json output
        # for unfiltered contigs

        unfiltered_total_length_value = 0
        unfiltered_number_of_contigs_value = 0
        unfiltered_average_covereage_value = 0
        unfiltered_N50_value = 0
        unfiltered_L50_value = 0
        unfiltered_number_of_Ns_value = 0

        # for filtered contigs

        filtered_total_length_value = 0
        filtered_number_of_contigs_value = 0
        filtered_average_covereage_value = 0
        filtered_N50_value = 0
        filtered_L50_value = 0
        filtered_number_of_Ns_value = 0

        with open(input_file_filtered) as f:
            for line in f:
                line = line.split("=")
                klucz, wartosc = line[0].rstrip(), line[1].rstrip()
                if klucz == "n_contig":
                    filtered_number_of_contigs_value = int(float(wartosc))
                elif klucz == "n_base":
                    filtered_total_length_value = int(float(wartosc))
                elif klucz == "ave_depth":
                    filtered_average_covereage_value = round(float(wartosc), 2)
                elif klucz == "n_N":
                    filtered_number_of_Ns_value = int(float(wartosc))
                elif klucz == "N50":
                    filtered_N50_value = int(float(wartosc))
                elif klucz == "L50":
                    filtered_L50_value = int(float(wartosc))

        # QC status podlega ocenia po uwzgledeniu przefiltrowanej listy contigow
        if filtered_N50_value < nvalue or \
            filtered_number_of_contigs_value > max_contigs or \
            filtered_average_covereage_value < coverage or \
            filtered_total_length_value < (float(completness) * int(genome_length)):
            status = "blad"
            QC_names = ['N50_value', "number_of_contigs_value", "average_coverage_value", "genome_completness"]
            QC_values = [filtered_N50_value,
                         filtered_number_of_contigs_value,
                         filtered_average_covereage_value,
                         filtered_total_length_value]

            QC_statuses = [filtered_N50_value < nvalue,
                           filtered_number_of_contigs_value > max_contigs,
                           filtered_average_covereage_value < coverage,
                           filtered_total_length_value < (float(completness) * int(genome_length))]
            error = f"Predicted genome does not meet criteria for following parameters:\t"
            for qc_status, qc_name, qc_value in zip(QC_statuses, QC_names ,QC_values):
                if qc_status:
                    error += f"{qc_name} is {qc_value}\t"
            filtered_json = {"step_name": "post-filtering",
                             "status": status,
                             "error_message": error}
        else:
            filtered_json = {"step_name": "post-filtering",
                             "status": "tak",
                             "total_length_value": filtered_total_length_value,
                             "number_of_contigs_value": filtered_number_of_contigs_value,
                             "average_coverage_value": filtered_average_covereage_value,
                             "N50_value": filtered_N50_value,
                             "L50_value": filtered_L50_value,
                             "number_of_Ns_value": filtered_number_of_Ns_value}

        with open(input_file_unfiltered) as f:
            for line in f:
                line = line.split("=")
                klucz, wartosc = line[0].rstrip(), line[1].rstrip()
                if klucz == "n_contig":
                    unfiltered_number_of_contigs_value = int(float(wartosc))
                elif klucz == "n_base":
                    unfiltered_total_length_value = int(float(wartosc))
                elif klucz == "ave_depth":
                    unfiltered_average_covereage_value = round(float(wartosc), 2)
                elif klucz == "n_N":
                    unfiltered_number_of_Ns_value = int(float(wartosc))
                elif klucz == "N50":
                    unfiltered_N50_value = int(float(wartosc))
                elif klucz == "L50":
                    unfiltered_L50_value = int(float(wartosc))

        unfiltered_json = {"step_name": "pre-filtering",
                           "status": "tak",
                           "total_length_value": unfiltered_total_length_value,
                           "number_of_contigs_value": unfiltered_number_of_contigs_value,
                           "average_coverage_value": unfiltered_average_covereage_value,
                           "N50_value": unfiltered_N50_value,
                           "L50_value": unfiltered_L50_value,
                           "number_of_Ns_value": unfiltered_number_of_Ns_value}
        json_output.append(filtered_json)
        json_output.append(unfiltered_json)

    with open(output, 'w') as f1:
        f1.write(json.dumps(json_output, indent = 4))
    print(status)
    return status


if __name__ == '__main__':
    # The main program returns 3 variables: "status", number_of_reads in a sample, and median_quality_of_reads
    # These can be used to determine QC status in tha module
    if len(sys.argv) == 1:
        main_program(['--help'])
    else:
        main_program(sys.argv[1:])
