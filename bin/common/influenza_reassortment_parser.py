#!/usr/bin/env python3

"""
This script parse the results of picard_wgsmetrics module in the viral pipeline to
produce a valid json (the viral_genom_data group of fileds)
"""

import sys
import click
import json
import re
import numpy as np

def parse_mapping_file(plik):
    """
    Skrypt do konwersji pliku z informacja jaki szczep i jego id jest referencja dla danego podtypu i kladu
    Plik ma postac macierzy gdzie 1sz wiersz to naglowek
    a Kolejne wiersze to nazwa "wewentrzna" genomu, a potem informacja jaka jest nazwa szczepu dla kazdego segmentu
    wynika to z faktu, ze dla jednego podtpyu H3N2_v nie dalo sie zlozyc genomu z jednego szczpeu
    @param plik: plik zmapowaniem
    @type plik: str
    @return: Slownik gdzie kluczem glownym jest nazwa genomu, podkluczami nazwy segmentu a wartosci dwuelemenotowa lista
    gdzie pierwszy element to nazwa szczepu a drugi jego identyfiaktor w gisaid
    @rtype: Dict[Dict -> List]
    """
    slownik_wynikow = {}
    with open(plik) as f:
        for line in f:
            line = line.split("\t")
            line[-1] = line[-1].rstrip()
            if line[0] == "#Directory":
                header = line[1:] # niech sam plik trzyma klucze podslownika (moze miznie na chrX_...)
            else:
                slownik_wynikow[line[0]] = {}

                # dane dla kazdego segmentu
                slownik_wynikow[line[0]][header[0]] = line[1].split(" ")
                slownik_wynikow[line[0]][header[1]] = line[2].split(" ")
                slownik_wynikow[line[0]][header[2]] = line[3].split(" ")
                slownik_wynikow[line[0]][header[3]] = line[4].split(" ")
                slownik_wynikow[line[0]][header[4]] = line[5].split(" ")
                slownik_wynikow[line[0]][header[5]] = line[6].split(" ")
                slownik_wynikow[line[0]][header[6]] = line[7].split(" ")
                slownik_wynikow[line[0]][header[7]] = line[8].split(" ")

    return slownik_wynikow

def parse_intermediate(plik, mapowania):
    segment_dict = {}
    list_of_refernces = []

    i = 0 # We only need to look at lines 0 and 1
    with open(plik) as f:
        for line in f:
            j = 0
            line = line.split()
            line[-1] = line[-1].rstrip()
            for segment in line:
                if i == 0:
                    segment_dict[str(j)] = [segment]
                if i == 1:
                    segment_dict[str(j)].append(segment)
                    segment_clean=segment.split('_')[0]

                    if re.findall('H5N\\d', segment_clean):
                        segment_clean="H5Nx" # We do not considere a reassortment between H5Nx as valid ones
                    list_of_refernces.append(segment_clean)
                j += 1
            i += 1


    reference_genome_data = []
    for element in segment_dict:
        segment_name, segment_reference = segment_dict[element]
        reference_genome_data.append({"segment_name": segment_name,
                                      "reference_subtype_name" : segment_reference.split('_')[0],
                                      "reference_strain_name" : f"{mapowania[segment_reference][segment_name][0]}",
                                      "reference_strain_id" : f"{mapowania[segment_reference][segment_name][1]}"})


    if len(np.unique(list_of_refernces)) > 1:
        reassortment_status="tak"
    else:
        reassortment_status = "nie"

    typ = ""
    for element in np.unique(list_of_refernces):
        if typ == "" and element not in ["Yamagata", "Victoria"]:
            typ="A"
        elif typ == "" and element in ["Yamagata", "Victoria"]:
            typ = "B"
        elif typ == "A" and element not in ["Yamagata", "Victoria"]:
            typ = "A"
        elif typ == "A" and element in ["Yamagata", "Victoria"]:
            typ = "unk"
        elif typ == "B" and element not in ["Yamagata", "Victoria"]:
            typ = "unk"
        elif typ == "B" and element in ["Yamagata", "Victoria"]:
            typ = "B"
        else:
            # this should result in an error when validating json schema
            typ="strange"

    return reference_genome_data, reassortment_status, typ


@click.command()
@click.option('-i', '--input_file', help='[INPUT] intermediate.txt file generated within a module',
              type=click.Path(), required=False)
@click.option('-u', '--subtype', help='[INPUT] Subtype predicted for this module',
              type=click.Path(), required=False)
@click.option('-s', '--status', help='[INPUT] PREDEFINED status that is transferred to an output json. '
                                     'If this status was either nie or blad fastqc will not run',
              type=click.Choice(['tak', 'nie', 'blad'], case_sensitive=False),  required=True)
@click.option('-m', '--mapping', help='[INPUT] PREDEFINED file with mapping between internal names of genomes'
                                      'like H3N2_3C2a1b2a2a3a to strains name like A/Denmark/65/2021 ',
              type=click.Path(), required=False)
@click.option('-r', '--error', help='[INPUT] PREDEFINED error message that is put in json. '
                                    'Only used when status was set to nie or blad',
              type=str,  required=False, default="")
@click.option('-o', '--output', help='[Output] Name of a file with json output',
              type=str,  required=True)

def main_program(status, output, input_file, mapping, subtype, error=""):
    if status != "tak":
        json_output = {"reference_genome_prediction_status": status,
                       "reference_genome_prediction_error_message": error}
    else:
        mapowania_genom_segment = parse_mapping_file(mapping)
        reference_genome_prediction_data, reasortment_result, typ = parse_intermediate(plik=input_file,
                                                                                       mapowania=mapowania_genom_segment)
        json_output = {
            "reference_genome_prediction_status": status,
            "reasortment_result": reasortment_result,
            "reference_genome_prediction_data": reference_genome_prediction_data,
            "type_name": typ,
            "subtype_name": subtype.split('_')[0]
        }
    with open(output, 'w') as f1:
        f1.write(json.dumps(json_output, indent=4))

    return True


if __name__ == '__main__':
    # The main program returns 3 variables: "status", number_of_reads in a sample, and median_quality_of_reads
    # These can be used to determine QC status in tha module
    if len(sys.argv) == 1:
        main_program(['--help'])
    else:
        main_program(sys.argv[1:])
