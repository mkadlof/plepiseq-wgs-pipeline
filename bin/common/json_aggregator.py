#!/usr/bin/env python3

import argparse
import datetime
import json


def json_aggregator(args):
    output = {"output": {}}
    output["output"]["pipeline_version"] = args.version
    output["output"]["pathogen"] = args.pathogen
    output["output"]["sampleId"] = args.sampleId
    output["output"]["created_timestamp"] = datetime.datetime.now().isoformat()

    output["output"]["dehumanized_fastq_data"] = {}
    output["output"]["viral_classification_data"] = {}
    output["output"]["viral_genome_data"] = {}
    output["output"]["viral_mutation_data"] = {}
    output["output"]["contamination_data"] = []
    output["output"]["genome_files_data"] = {}
    output["output"]["sequencing_summary_data"] = []

    if output["output"]["pathogen"] == "sars2":
        output["output"]["sars_data"] = {}

    with open("output.json", "w") as f:
        json.dump(output, f, indent=4)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('version', help="Pipeline version")
    parser.add_argument('pathogen', help="Pathogen")
    parser.add_argument('sampleId', help="Sample ID")
    parser.add_argument('--wgsMetrics', help="WGS metrics file")
    args = parser.parse_args()

    json_aggregator(args)


if __name__ == '__main__':
    main()
