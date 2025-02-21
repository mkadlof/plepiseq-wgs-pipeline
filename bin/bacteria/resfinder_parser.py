
import sys
import click
import json


@click.command()
@click.option('-i', '--input_file_resfinder', help='[INPUT] a path to an input file with Resfinder results',
              type=click.Path(), default="ResFinder_results_tab.txt", required=True)
@click.option('-j', '--input_file_pointfinder', help='[INPUT] a path to an input file with Pointfinder results',
              type=click.Path(), default="PointFinder_results.txt", required=True)
@click.option('-s', '--status', help='[INPUT] PREDEFINED status that is transferred to an output json. '
                                     'If this status was either nie or blad fastqc will not run',
              type=click.Choice(['tak', 'nie', 'blad'], case_sensitive=False),  required=True)
@click.option('-r', '--error', help='[INPUT] PREDEFINED error message that is put in json. '
                                    'Only used when status was set to nie or blad',
              type=str,  required=False, default="")
@click.option('-o', '--output', help='[Output] Name of a file with json output',
              type=str,  required=True)
def main_program(status, input_file_resfinder, input_file_pointfinder, output, error=""):
    if status != "tak":
        json_output = {"program_name" : "ResFinder/PointFinder",
                        "status": status,
                        "error_message": error}
    else:
        json_output = {"program_name" : "ResFinder/PointFinder",
                     "status": status}
        slownik = {}
        with open(input_file_resfinder) as f:
            for line in f:
                line = line.split("\t")
                line[-1] = line[-1].rstrip()
                if "Resistance gene" == line[0]:
                    continue
                else:
                    resistance_gene, seq_identity, coverage, contig_name, antibiotic, reference_name = line[0], line[1], line[3], line[5], line[7], line[8]
                    antibiotic_list = antibiotic.split(",")
                    for element in antibiotic_list:
                        if element not in slownik.keys():
                            slownik[element] = []
                        slownik[element].append([resistance_gene, contig_name, seq_identity, coverage, reference_name, "gen"])
        with open(input_file_pointfinder) as f:
            for line in f:
                line = line.split("\t")
                line[-1] = line[-1].rstrip()
                if "Mutation" == line[0]:
                    continue
                else:
                    mutation,  antibiotic = line[0], line[3]
                    antibiotic_list = antibiotic.split(",")
                    gene_name, mutation_name = mutation.split(" ")
                    for element in antibiotic_list:
                        if element not in slownik.keys():
                            slownik[element] = []
                        slownik[element].append([gene_name, mutation_name, "mutacja_punktowa"])
        # skladanie wlasciwego jsona
        json_output["program_data"] = []
        for antibiotic, czynniki in slownik.items():
            tmp_dict = {}
            tmp_dict["antibiotic_name"] = antibiotic
            tmp_dict["antibiotic_status"] = "oporny"
            tmp_dict["antibiotic_resistance_data"] = []
            for czynnik in czynniki:
                if czynnik[-1] == "gen":
                    tmp_dict["antibiotic_resistance_data"].append(
                        {"factor_name" : czynnik[0],
                         "factor_contig_name" : czynnik[1],
                         "factor_sequence_similarity_to_reference_value": int(float(czynnik[2])),
                         "factor_degree_of_overlap_with_reference_value" : int(float(czynnik[3])),
                         "factor_reference_name" : czynnik[4],
                         "factor_type_name" : "gen"
                         }
                    )
                elif czynnik[-1] == "mutacja_punktowa":
                    tmp_dict["antibiotic_resistance_data"].append(
                        {"factor_name": czynnik[0],
                         "factor_mutation": czynnik[1],
                         "factor_type_name": "mutacja_punktowa"
                         }
                    )
            json_output["program_data"].append(tmp_dict)

    with open(output, 'w') as f1:
        f1.write(json.dumps(json_output, indent = 4))

    return True

if __name__ == '__main__':
    # The main program returns 3 variables: "status", number_of_reads in a sample, and median_quality_of_reads
    # These can be used to determine QC status in tha module
    if len(sys.argv) == 1:
        main_program(['--help'])
    else:
        main_program(sys.argv[1:])

