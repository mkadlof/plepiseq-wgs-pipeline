
import sys
import click
import json


@click.command()
@click.option('-i', '--input_file', help='[INPUT] a path to an input file with VFDB results',
              type=click.Path(),  required=True)
@click.option('-s', '--status', help='[INPUT] PREDEFINED status that is transferred to an output json. '
                                     'If this status was either nie or blad fastqc will not run',
              type=click.Choice(['tak', 'nie', 'blad'], case_sensitive=False),  required=True)
@click.option('-r', '--error', help='[INPUT] PREDEFINED error message that is put in json. '
                                    'Only used when status was set to nie or blad',
              type=str,  required=False, default="")
@click.option('-o', '--output', help='[Output] Name of a file with json output',
              type=str,  required=True)
def main_program(status, input_file, output, error=""):
    if status != "tak":
        json_output = {"status": status,
                       "program_name" : "VFDB",
                        "error_message": error}
    else:
        json_output = {"status": status,
                       "program_name" : "VFDB",
                       "patotype_data" : [{"patotype_name" : "NA"}]}

        with open(input_file) as f:
            virulence_genes_data = []
            for line in f:
                line = line.split("\t")
                line[-1] = line[-1].rstrip()
                gene_name, seq_id, coverage, contig_name, reference_accession_name = line[3], line[5], line[6], line[4], line[7]
                if line[-1] == "BRAK":
                    pass
                else:
                    virulence_genes_data.append({"gene_name" : gene_name,
                                                 "sequence_similarity_to_reference_value" : round(float(seq_id), 2),
                                                 "degree_of_overlap_with_reference_value" : round(float(coverage), 2),
                                                 "contig_name" : contig_name,
                                                 "reference_accession_name" : reference_accession_name})
        json_output["virulence_genes_data"] = virulence_genes_data

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

