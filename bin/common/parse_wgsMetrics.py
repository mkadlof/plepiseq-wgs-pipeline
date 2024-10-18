#!/usr/bin/env python3

import argparse
import json


def parse_file(input: str, output: str):
    with open(input, 'r') as file:
        lines = file.readlines()
        mean_coverage = float(lines[7].split()[1])

    with open(output, 'w') as file:
        json.dump({'average_coverage_value': mean_coverage}, file, indent=4)


def main():
    parser = argparse.ArgumentParser(description='Parse wgsMetrics output')
    parser.add_argument('picard_stats', help='Input file (picard stats)')
    parser.add_argument('output', help='Output file')
    args = parser.parse_args()

    parse_file(args.picard_stats, args.output)


if __name__ == '__main__':
    main()
