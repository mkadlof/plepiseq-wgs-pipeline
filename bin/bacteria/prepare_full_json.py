import sys
import click
import json
import datetime

@click.command()
@click.option('-q', '--sistr_file', help='[INPUT] a path to an input file with sistr json',
              type=click.Path(),  required=True)
@click.option('-w', '--seqsero_file', help='[INPUT] a path to an input file with seqsero json',
              type=click.Path(),  required=True)
@click.option('-e', '--spifinder_file', help='[INPUT] a path to an input file with spifinder json',
              type=click.Path(),  required=True)
@click.option('-r', '--ectyper_file', help='[INPUT] a path to an input file with ectyper json',
              type=click.Path(),  required=True)
@click.option('-t', '--virulencefinder_file', help='[INPUT] a path to an input file with virulencefinder json',
              type=click.Path(),  required=True)
@click.option('-y', '--vfdb_file', help='[INPUT] a path to an input file with VFDB json',
              type=click.Path(),  required=True)
@click.option('-u', '--plasmidfinder_file', help='[INPUT] a path to an input file with plasmidfinder json',
              type=click.Path(),  required=True)
@click.option('-i', '--amrfinder_file', help='[INPUT] a path to an input file with amrfinder json',
              type=click.Path(),  required=True)
@click.option('-o', '--resfinder_file', help='[INPUT] a path to an input file with resfinder json',
              type=click.Path(),  required=True)
@click.option('-p', '--cgmlst_file', help='[INPUT] a path to an input file with cgMLST json',
              type=click.Path(),  required=True)
@click.option('-a', '--mlst_file', help='[INPUT] a path to an input file with MLST json',
              type=click.Path(),  required=True)
@click.option('-s', '--fastqc_forward_file', help='[INPUT] a path to an input file with json for'
                                                  ' forward fastq file',
              type=click.Path(),  required=True)
@click.option('-d', '--fastqc_reverse_file', help='[INPUT] a STRING with a PATH to an input file with '
                                                  'json for reverse fastq file. Can be "skip". In that case'
                                                  ' this parameter is not used to build final json' ,
              type=str,  required=False)
@click.option('-f', '--contaminations_file', help='[INPUT] a path to an input file with json containing'
                                                  ' summary of get_specied workflow',
              type=click.Path(),  required=True)
@click.option('--genus_species_file', help='[INPUT] a path to an input file with json containing'
                                                  'infor regading genus and species predicted for a sample',
              type=click.Path(),  required=True)
@click.option('-g', '--initial_mlst_file', help='[INPUT] a path to an input file with json for mlst '
                                                'prediction on raw fastq files',
              type=click.Path(),  required=True)
@click.option('-j', '--genome_statistics_file', help='[INPUT] a path to an input file with json with '
                                                     'statistics for predicted genome',
              type=click.Path(),  required=True)
@click.option('-k', '--genome_file', help='[INPUT] a path to an input file with json with paths to fastas '
                                          'for each segment',
              type=click.Path(),  required=True)
@click.option('-m', '--patotyp', help='[INPUT] String with all patotypes names predicted for a sample',
              type=str,  required=True)
@click.option('--alphafold_file', help='[INPUT] Json prodyced by the alphafold module',
              type=click.Path(),  required=True)
@click.option('-w', '--executiondir', help='[INPUT] String with execution path of the program',
              type=str,  required=True)
@click.option('-l', '--repo_version', help='[INPUT] Version of a repository',
              type=str,  required=True)
@click.option('-z', '--output', help='[Output] Name of a file with json output',
              type=str,  required=True)
def main_program(sistr_file, seqsero_file, spifinder_file,  ectyper_file,  virulencefinder_file,  alphafold_file, vfdb_file,
                 plasmidfinder_file, amrfinder_file, resfinder_file,  cgmlst_file, mlst_file,  fastqc_forward_file,
                 fastqc_reverse_file, contaminations_file, genus_species_file, initial_mlst_file, genome_file, patotyp,
                 genome_statistics_file, repo_version, executiondir, output):

    json_output = {}
    json_output["output"] = {}
    json_output["output"]["pathogen"] = "bacteria"
    json_output["output"]["pathogen_predicted_genus"] = json.load(open(genus_species_file))['pathogen_predicted_genus']
    json_output["output"]["pathogen_predicted_species"] = json.load(open(genus_species_file))['pathogen_predicted_species']
    json_output["output"]["ExecutionDir_dir"] = executiondir
    json_output["output"]["timestamp"] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    json_output["output"]["pipeline_version"] = repo_version
    json_output["output"]["genome_files_data"] = json.load(open(genome_file))
    if fastqc_reverse_file != "skip":
        json_output["output"]["sequencing_summary_data"] = [json.load(open(fastqc_forward_file))[0],
                                                       json.load(open(fastqc_reverse_file))[0]]
    else:
        json_output["output"]["sequencing_summary_data"] = [json.load(open(fastqc_forward_file))[0]]

    json_output["output"]["contamination_data"] = json.load(open(contaminations_file))
    json_output["output"]["initial_mlst_data"] = json.load(open(initial_mlst_file))
    json_output["output"]["bacterial_genome_data"] = json.load(open(genome_statistics_file))
    json_output["output"]["antigenic_data"] = [json.load(open(sistr_file)),
                                          json.load(open(seqsero_file)),
                                          json.load(open(ectyper_file))]
    json_output["output"]["virulence_islands_data"] = json.load(open(spifinder_file))
    json_output["output"]["structural_data"] = json.load(open(alphafold_file))

    json_output["output"]["virulence_genes_data"] = [json.load(open(virulencefinder_file))]

    # Dla VFDB musimy jeszcze dodac wyniki modulu parse_VFDB_ecoli, ktory dla ecoli okresla poprawny patotyp
    # Dla innych organizmow daje patotyp NA
    vfdb_output = json.load(open(vfdb_file))
    if vfdb_output["status"] == "tak":
        del(vfdb_output["patotype_data"]) # wartosc tworzona aby mozna bylo poprawnie ealuwac czastkowy schemat
        # zawieral dummy wartosc
        vfdb_output["patotype_data"] = []
        patotyp = patotyp.split(" ") # Patotyp to string ktory po spacji zawierac moze wiele patotypow np. "STEC EIEC"
        for pato in patotyp:
            vfdb_output["patotype_data"].append({"patotype_name": pato})

    json_output["output"]["virulence_genes_data"].append(vfdb_output)
    json_output["output"]["plasmids_data"] = json.load(open(plasmidfinder_file))
    json_output["output"]["drug_resistance_data"] = [json.load(open(amrfinder_file)),
                                                json.load(open(resfinder_file))]
    json_output["output"]["mlst_data"] = [json.load(open(mlst_file)),
                                     json.load(open(cgmlst_file))]

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
