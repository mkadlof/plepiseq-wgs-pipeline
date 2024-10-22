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
            # filename example: segment_MN908947.3.bedgraph
            segment_name = line.split("_")[1].split(".")[0]
            data_file_path = args.publish_dir + "/" + line.strip()
            output["output"]["viral_genome_data"]["coverage_barplot_data"].append(
                {"segment_name": segment_name, "data_file_path": data_file_path})
    if args.pathogen == "sars2":
        output["output"]["viral_genome_data"][
            "coverage_histogram_data"] = [{
            "segment_name": "MN908947.3",
            "data_file_path": args.publish_dir + "/" + "coverage_histogram_data.tsv"
        }]
    else:
        if args.pathogen == "sars2":
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
    if output["output"]["pathogen"] == "sars2":
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


def json_aggregator(args):
    output = {"output": {}}
    output["output"]["pipeline_version"] = args.version
    output["output"]["pathogen"] = args.pathogen
    output["output"]["sampleId"] = args.sampleId
    output["output"]["timestamp"] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    output["output"]["dehumanized_fastq_data"] = []
    output["output"]["viral_classification_data"] = []
    output["output"]["viral_genome_data"] = {}
    output["output"]["viral_mutation_data"] = []
    output["output"]["contamination_data"] = []
    output["output"]["genome_files_data"] = {}
    output["output"]["sequencing_summary_data"] = []

    # fill json sections
    fill_viral_genome_data(args, output)
    fill_sars_data(args, output)

    output["output"]["genome_files_data"]['file_data'] = []  # TODO
    output["output"]["genome_files_data"]['status'] = "tak"  # TODO

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
    args = parser.parse_args()

    json_aggregator(args)


if __name__ == '__main__':
    main()
