#!/usr/bin/env python3

import argparse
import datetime
import json

from warnings import warn


def fill_viral_genome_data(args, output):
    # viral_genome_data
    with open(args.wgsMetrics) as f:
        wgsMetrics = json.load(f)
        average_coverage_value = round(wgsMetrics["average_coverage_value"], 2)
        output["output"]["viral_genome_data"]["average_coverage_value"] = average_coverage_value
    output["output"]["viral_genome_data"]["coverage_barplot_data"] = []
    with open(args.segment_bedgraphs_files) as f:
        lines = f.readlines()
        for line in lines:
            # filename example: segment_MN908947.3.bedgraph, segment_chr4_HA.bedgraph
            segment_name = line.split("_")[-1].split(".")[0]
            data_file_path = args.publish_dir + "/" + line.strip()
            output["output"]["viral_genome_data"]["coverage_barplot_data"].append(
                {"segment_name": segment_name, "data_file_path": data_file_path})
    if args.pathogen == "sars2":
        output["output"]["viral_genome_data"][
            "coverage_histogram_data"] = [{
            "segment_name": "MN908947.3",
            "data_file_path": args.publish_dir + "/" + "coverage_histogram_data.tsv"
        }]
    elif args.pathogen == "influenza":
        output["output"]["viral_genome_data"][
            "coverage_histogram_data"] = [{
            "segment_name": "",
            "data_file_path": ""  # TODO:
        }]
        warn(
            "coverage_histogram_data is not implemented for pathogen: " + args.pathogen + ". Fields in output.json are empty.")
    with open(args.consensus) as f:
        consensus = json.load(f)
        output["output"]["viral_genome_data"]['total_length_value'] = consensus['total_length_value']
        output["output"]["viral_genome_data"]['number_of_Ns_value'] = consensus['number_of_Ns_value']
    output["output"]["viral_genome_data"]["primer_usage_data"] = []  # TODO
    output["output"]["viral_genome_data"]["error_message"] = ""  # TODO
    output["output"]["viral_genome_data"]["status"] = "tak"  # TODO
    # end viral_genome_data


def fill_sars_data(args, output):
    output["output"]["sars_data"] = {
        'coinfection_histogram_file': "",  # TODO
        'coinfection_message': "",  # TODO
        'coinfections_status': "nie",  # TODO
        'freyja_lineage1_abundance': -1,  # TODO
        'freyja_lineage1_name': "",  # TODO
        'freyja_lineage2_abundance': -1,  # TODO
        'freyja_lineage2_name': "",  # TODO
        'protein_structure_data': []  # TODO
    }


def fill_infl_data(args, output):
    output["output"]["infl_data"] = {
        'resistance_data': [],  # TODO
        'protein_structure_data': [],  # TODO
        'reassortment_data': [],  # TODO
        'subtype_name': "...",  # TODO
        'type_name': "unk",  # TODO
    }


def fill_contamination_data(args, output):
    with open(args.contamination) as f:
        contamination = json.load(f)
        output["output"]["contamination_data"] = contamination


def fill_sequencing_summary_data(args, output):
    files = [args.fastqc_pre[0], args.fastqc_pre[1], args.fastqc_post[0], args.fastqc_post[1]]
    for file in files:
        with open(file) as f:
            fastqc = json.load(f)
            output["output"]["sequencing_summary_data"].append(fastqc[0])


def fill_genome_files_data(args, output):
    output["output"]['genome_files_data']['file_data'] = []
    output["output"]["genome_files_data"]['status'] = "tak"  # TODO
    if args.pathogen == 'sars2':
        # For sars file is hardcoded. It never changes.
        output["output"]['genome_files_data']['file_data'].append({
            'segment_name': 'MN908947.3',
            'segment_path': args.publish_dir + '/output_consensus_masked_SV.fa'
        })
    elif args.pathogen == 'influenza':
        # For influenza, we need to read the list of fasta files from a file
        with open(args.list_of_fasta_files) as f:
            for line in f:
                filename = line.strip()
                segment_name = filename.split("_")[1].split(".")[0]
                output["output"]['genome_files_data']['file_data'].append({
                    'segment_name': segment_name,
                    'segment_path': filename
                })


def fill_viral_classification_data(file_path, output):
    if file_path == "non-existent":
        return
    with open(file_path) as f:
        data = json.load(f)
        print(data)
        for i in data:
            output["output"]["viral_classification_data"].append(i)


def json_aggregator(args):
    output = {"output": {}}

    # Fields independent of modules output
    output["output"]["pipeline_version"] = args.version
    output["output"]["pathogen"] = args.pathogen
    output["output"]["sampleId"] = args.sampleId
    output["output"]["timestamp"] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    # Initialize empty fields
    output["output"]["dehumanized_fastq_data"] = []
    output["output"]["viral_classification_data"] = []
    output["output"]["viral_genome_data"] = {}
    output["output"]["viral_mutation_data"] = []
    output["output"]["contamination_data"] = []
    output["output"]["genome_files_data"] = {}
    output["output"]["sequencing_summary_data"] = []

    # fill json dehumanized sections
    output["output"]["dehumanized_fastq_data"] = [
        args.publish_dir + "/forward_paired_nohuman.fastq.gz",
        args.publish_dir + "/reverse_paired_nohuman.fastq.gz",
    ]

    # fill contamination data (from kraken2 module)
    fill_contamination_data(args, output)

    # paths to fasta files
    fill_genome_files_data(args, output)

    # fill json sections that depend on pathogen
    fill_viral_genome_data(args, output)
    if output["output"]["pathogen"].lower() in ["sars2", "sars-cov-2"]:
        fill_sars_data(args, output)
    elif output["output"]["pathogen"].lower() in ["influenza"]:
        fill_infl_data(args, output)
    elif output["output"]["pathogen"].lower() in ["rsv"]:
        raise NotImplementedError("RSV is not implemented yet")

    fill_sequencing_summary_data(args, output)

    fill_viral_classification_data(args.pangolin, output)
    fill_viral_classification_data(args.nextclade, output)

    with open("output.json", "w") as f:
        json.dump(output, f, indent=4)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('version', help="Pipeline version")
    parser.add_argument('pathogen', help="Pathogen")
    parser.add_argument('sampleId', help="Sample ID")
    parser.add_argument('publish_dir', help="Publish directory")
    parser.add_argument('--wgsMetrics', help="WGS metrics file")
    parser.add_argument('--segment_bedgraphs_files', help="Segment bedgraphs files")
    parser.add_argument('--consensus', help="JSON from consensus module")
    parser.add_argument('--list_of_fasta_files', help="File with a list of fastq files")
    parser.add_argument('--contamination', help="JSON from contamination detection (kraken2)")
    parser.add_argument('--fastqc_pre', nargs=2, help='Two json files from fastqc module (pre_filtering)')
    parser.add_argument('--fastqc_post', nargs=2, help='Two json files from fastqc module (post_filtering)')
    parser.add_argument('--pangolin', help="JSON from viral classification module (pangolin)")
    parser.add_argument('--nextclade', help="JSON from viral classification module (nextclade)")

    args = parser.parse_args()

    json_aggregator(args)


if __name__ == '__main__':
    main()
