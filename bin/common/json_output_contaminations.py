#!/usr/bin/env python3
import json
import sys
import click


@click.command()
@click.option('-k', '--input_kraken', help="[INPUT] a path to a kraken2 output. Use "
                                           "\"skip\" to not include results of this program in json output",
              type=str,  required=True)
@click.option('-g', '--input_metaphlan_genera', help="[INPUT] a path to a metaphlan output. Use "
                                                     "\"skip\" to not include results of this program in json output",
              type=str,  required=True)
@click.option('-x', '--input_metaphlan_species', help="[INPUT] a path to a metaphlan output. Use "
                                                      "\"skip\" to not include results of this program in json output",
              type=str,  required=True)
@click.option('-y', '--input_kmerfinder', help="[INPUT] a path to a kmerfinder output. Use"
                                               "\"skip\" to not include results of this program in json output ",
              type=str,  required=True)
@click.option('-s', '--status', help='[INPUT] PREDEFINED status that is transferred to an output json. '
                                     'If this status was either nie or blad fastqc will not run',
              type=click.Choice(['tak', 'nie', 'blad'], case_sensitive=False),  required=True)
@click.option('-m', '--error', help='[INPUT] PREDEFINED error message that is put in json. '
                                    'Only used when status was set to nie or blad',
              type=str,  required=False, default="")
@click.option('-o', '--output', help='[Output] Name of a file with json output',
              type=str,  required=True)
def main_program(input_kraken, input_metaphlan_genera, input_metaphlan_species, input_kmerfinder,
                 status, error, output):
    full_output = []
    # Kraken2 results section
    if input_kraken == "skip":
        # ignore kraken2 whatsoever
        pass
    else:
        if status != "tak":
            kraken2_json = {"program_name": "kraken2",
                            "status": status,
                            "error_message": error}
        else:
            f = open(input_kraken).readlines()
            genus_dict = {}
            species_dict = {}
            for line in f:
                line = line.split()
                if line[3] == "S":
                    try:
                        species_dict[" ".join(list(map(str, line[5:])))] = float(line[0])
                        #species_dict[line[5] + " " + line[6]] = float(line[0])
                    except IndexError:
                        species_dict[line[5]] = float(line[0])
                elif line[3] == "G":
                    genus_dict[line[5]] = float(line[0])
            genus_names_sorted = sorted(genus_dict, key=lambda x: genus_dict[x], reverse=True)
            species_names_sorted = sorted(species_dict, key=lambda x: species_dict[x], reverse=True)
            kraken2_json = {"program_name": "kraken2",
                            "status": "tak",
                            "main_genus_name": genus_names_sorted[0],
                            "secondary_genus_name": genus_names_sorted[1],
                            "main_species_name": species_names_sorted[0],
                            "secondary_species_name": species_names_sorted[1],
                            "main_genus_value": round(genus_dict[genus_names_sorted[0]], 2),
                            "secondary_genus_value": round(genus_dict[genus_names_sorted[1]], 2),
                            "main_species_value": round(species_dict[species_names_sorted[0]], 2),
                            "secondary_species_value": round(species_dict[species_names_sorted[1]], 2)}

        full_output.append(kraken2_json)

    # Metaphlan results section
    if input_metaphlan_genera == "skip" and input_metaphlan_species == "skip":
        pass
    else:
        if status != "tak":
            metaphlan_json = {"program_name": "metaphlan",
                              "status": status,
                              "error_message": error}
        else:
            f1 = open(input_metaphlan_genera).readlines()
            f2 = open(input_metaphlan_species).readlines()
            genus_dict = {}
            species_dict = {}
            for line in f1:
                line = line.split()
                if "g__" in line[0]:
                    genus_dict[line[0].split('_')[-1]] = float(line[2])
            for line in f2:
                line = line.split()
                if "s__" in line[0]:
                    species_dict[" ".join(line[0].split('_')[-2:])] = float(line[2])

            genus_names_sorted = sorted(genus_dict, key=lambda x: genus_dict[x], reverse=True)
            species_names_sorted = sorted(species_dict, key=lambda x: species_dict[x], reverse=True)
            if len(genus_names_sorted) == 1:
                genus_names_sorted.append('None')
                genus_dict['None'] = 0
            if len(species_names_sorted) == 1:
                species_names_sorted.append('None')
                species_dict['None'] = 0

            metaphlan_json = {"program_name": "metaphlan",
                              "status": "tak",
                              "main_genus_name": genus_names_sorted[0],
                              "secondary_genus_name": genus_names_sorted[1],
                              "main_species_name": species_names_sorted[0],
                              "secondary_species_name": species_names_sorted[1],
                              "main_genus_value": round(genus_dict[genus_names_sorted[0]], 2),
                              "secondary_genus_value": round(genus_dict[genus_names_sorted[1]], 2),
                              "main_species_value": round(species_dict[species_names_sorted[0]], 2),
                              "secondary_species_value": round(species_dict[species_names_sorted[1]], 2)}

        full_output.append(metaphlan_json)
    # Kmerfinder
    if input_kmerfinder == "skip":
        pass
    else:
        if status != "tak":
            kmerfinder_json = {"program_name": "kmerfinder",
                               "status": status,
                               "error_message": error}
        else:
            f = open(input_kmerfinder).readlines()
            genus_dict = {}
            species_dict = {}
            for line in f:
                line = line.split("\t")
                if "#" in line[0]:
                    continue
                rodzaj = line[14].split(' ')[0]
                gatunek = rodzaj + " " + line[14].split(' ')[1]
                pokrycie = float(line[8])

                if rodzaj in genus_dict.keys():
                    genus_dict[rodzaj] = genus_dict[rodzaj] + pokrycie
                else:
                    genus_dict[rodzaj] = pokrycie

                if gatunek in genus_dict.keys():
                    species_dict[gatunek] = species_dict[gatunek] + pokrycie
                else:
                    species_dict[gatunek] = pokrycie

            genus_names_sorted = sorted(genus_dict, key=lambda x: genus_dict[x], reverse=True)
            species_names_sorted = sorted(species_dict, key=lambda x: species_dict[x], reverse=True)
            if len(genus_names_sorted) == 1:
                genus_names_sorted.append('None')
                genus_dict['None'] = 0
            if len(species_names_sorted) == 1:
                species_names_sorted.append('None')
                species_dict['None'] = 0

            kmerfinder_json = {"program_name": "kmerfinder",
                               "status": "tak",
                               "main_species_name": species_names_sorted[0],
                               "secondary_species_name": species_names_sorted[1],
                               "main_species_coverage": round(species_dict[species_names_sorted[0]], 2),
                               "secondary_species_coverage": round(species_dict[species_names_sorted[1]], 2)}


        full_output.append(kmerfinder_json)

    with open(output, 'w') as f:
        f.write(json.dumps(full_output, indent = 4))

    return True


if __name__ == '__main__':
    if len(sys.argv) == 1:
        main_program(['--help'])
    else:
        print(main_program(sys.argv[1:]))
