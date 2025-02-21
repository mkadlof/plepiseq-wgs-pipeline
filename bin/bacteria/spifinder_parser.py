
import sys
import click
import json


@click.command()
@click.option('-i', '--input_file', help='[INPUT] a path to an input file with spifinder results',
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
                        "error_message": error}
    else:
        json_output = {"status": status}
        slownik = {}
        with open(input_file) as f:
            spi_data = []
            for line in f:
                line = line.split("\t")
                line[-1] = line[-1].rstrip()
                if "Database" == line[0]:
                    continue
                else:
                    spi_name, seq_identity, coverage, contig_name, reference_name = line[1], int(float(line[2])), line[3], line[4], line[-1]

                    coverage_template, coverage_query = coverage.split("/")
                    coverage = float(coverage_template)/ float(coverage_query) * 100
                    if coverage > 100:
                        #  w output spifinder nie mamy dostepy do granic mapowan, tylko info jaka jest dlugosc query
                        #  stad jesli query ma insercje wzgeledm template jego dlugosc moze byc wieksza od dlugosci
                        #  template i dostaniemy wartosc > 100. Json tego nie zaakceptuje wiec w takiej sytuacji na
                        # na sztywno daje 100
                        coverage = 100
                    spi_data.append({"spi_name" : spi_name,
                                      "sequence_similarity_to_reference_value" : seq_identity,
                                      "degree_of_overlap_with_reference_value" : round(coverage, 2),
                                      "contig_name" : contig_name,
                                      "reference_accession_name" : reference_name})
        json_output["spi_data"] = spi_data

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

