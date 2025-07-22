#!/usr/bin/env python3

"""
This script parse the results of picard_wgsmetrics module in the viral pipeline to
produce a valid json (the viral_genom_data group of fileds)
"""

import sys
import click
import json


def parse_picard(plik, coverage_histogram_out):
    average_coverage = 0
    valid_line_coverage = False
    valid_line_histogram = False
    i = 1
    with open(plik) as f, open(coverage_histogram_out, 'w') as f1:
        for line in f:
            if "GENOME_TERRITORY" in line.split():
                valid_line_coverage = True
                continue
            if valid_line_coverage:
                average_coverage = line.split()[1]
                valid_line_coverage = False

            if "coverage" in line.split():
                valid_line_histogram = True
                f1.write(f'#indeks,pokrycie,licznosc\n')
                continue

            if valid_line_histogram:
                try:
                    line = line.split()
                    f1.write(f'{i},{line[0]},{line[1]}\n')
                    i += 1
                except IndexError:
                    valid_line_histogram = False

    return float(average_coverage)


@click.command()
@click.option('-i', '--input_file_picard', help='[INPUT] a path to the results of picatd wgsmetrics command',
              type=click.Path(), required=False)
@click.option('-j', '--input_file_primers', help='[INPUT] a path to Primer_usage.txt file produced by'
                                                 ' my filtering script',
              type=click.Path(), required=False)
@click.option('-l', '--input_file_bedgraph', help='[INPUT] a file with a list of files used for plotting '
                                                  'coverage BARPLOTS',
              type=click.Path(), required=False)
@click.option('-k', '--output_path', help='[INPUT] a path where all the relevant files',
              type=str)
@click.option('-s', '--status', help='[INPUT] PREDEFINED status that is transferred to an output json. '
                                     'If this status was either nie or blad fastqc will not run',
              type=click.Choice(['tak', 'nie', 'blad'], case_sensitive=False),  required=True)
@click.option('-r', '--error', help='[INPUT] PREDEFINED error message that is put in json. '
                                    'Only used when status was set to nie or blad',
              type=str,  required=False, default="")
@click.option('-o', '--output', help='[Output] Name of a file with json output',
              type=str,  required=True)
def main_program(status, output, input_file_picard, input_file_primers, input_file_bedgraph, output_path, error=""):
    if status != "tak":
        json_output = {"status": status,
                       "error_message": error}
    else:
        coverage_histogram_file = "coverage_histogram.csv"
        average_coverage = parse_picard(plik=input_file_picard,
                                        coverage_histogram_out=coverage_histogram_file)
        coverage_barplot_data = []
        with open(input_file_bedgraph) as f:
            for line in f:
                segment_name, segment_file = line.split(',')
                coverage_barplot_data.append({"segment_name": segment_name,
                                              "segment_file": f'{output_path}/{segment_file.rstrip()}'})

        primer_usage_data = []

        with open(input_file_primers) as f:
            for line in f:
                if "Segment" in line.split():
                    slownik_segmentu = {"segment_name": "dummy"}
                    continue

                segment_name, primer_name, primer_usage = line.rsplit()

                if segment_name != slownik_segmentu["segment_name"]:
                    # nowy segment w primerach "zzrzucamy stary dict do wynikow
                    if slownik_segmentu["segment_name"] != "dummy":
                        primer_usage_data.append(slownik_segmentu)
                    # zakladamy nowy pusty slownik
                    slownik_segmentu = {"segment_name": segment_name, "primer_data": []}

                slownik_segmentu["primer_data"].append({"primer_name": primer_name,
                                                        "primer_usage": int(primer_usage.rstrip())})

            # Zrzucamy tez ostatni segment
            primer_usage_data.append(slownik_segmentu)

        json_output = {
            "status": status,
            "average_coverage_value": float(f'{average_coverage:.2f}'),
            "coverage_histogram_file": f'{output_path}/{coverage_histogram_file}',
            "coverage_barplot_data": coverage_barplot_data,
            "primer_usage_data": primer_usage_data
        }
    with open(output, 'w') as f1:
        f1.write(json.dumps(json_output, indent=4, ensure_ascii=False))

    return True


if __name__ == '__main__':
    # The main program returns 3 variables: "status", number_of_reads in a sample, and median_quality_of_reads
    # These can be used to determine QC status in tha module
    if len(sys.argv) == 1:
        main_program(['--help'])
    else:
        main_program(sys.argv[1:])
