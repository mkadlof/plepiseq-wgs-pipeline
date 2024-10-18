#!/usr/bin/env python3

import argparse
import datetime
import json


def json_aggregator(args):
    output = {"output": {}}
    output["output"]["pipeline_version"] = args.version
    output["output"]["pathogen"] = args.pathogen
    output["output"]["sampleId"] = args.sampleId
    output["output"]["timestamp"] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    output["output"]["dehumanized_fastq_data"] = {}
    output["output"]["viral_classification_data"] = {}
    output["output"]["viral_genome_data"] = {}
    output["output"]["viral_mutation_data"] = {}
    output["output"]["contamination_data"] = []
    output["output"]["genome_files_data"] = {}
    output["output"]["sequencing_summary_data"] = []

    if output["output"]["pathogen"] == "sars2":
        output["output"]["sars_data"] = {}

    # viral_genome_data
    with open(args.wgsMetrics) as f:
        wgsMetrics = json.load(f)
        average_coverage_value = wgsMetrics["average_coverage_value"]
        output["output"]["viral_genome_data"]["average_coverage_value"] = average_coverage_value

    output["output"]["viral_genome_data"]["coverage_barplot_data"] = []
    with open(args.segment_bedgraphs_files) as f:
        lines = f.readlines()
        for line in lines:
            # filename example: segment_MN908947.3.bedgraph
            segment_name = line.split("_")[1].split(".")[0]
            data_file_path = args.publish_dir + "/" + line.strip()
            output["output"]["viral_genome_data"]["coverage_barplot_data"].append({"segment_name": segment_name, "data_file_path": data_file_path})
    output["output"]["viral_genome_data"]["coverage_histogram_data"] = args.publish_dir + "/" + "coverage_histogram_data.tsv"

    with open(args.consensus) as f:
        consensus = json.load(f)
        output["output"]["viral_genome_data"]['total_length_value'] = consensus['total_length_value']
        output["output"]["viral_genome_data"]['number_of_Ns_value'] = consensus['number_of_Ns_value']

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
