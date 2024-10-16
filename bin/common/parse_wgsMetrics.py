#!/usr/bin/env python3

import argparse
import json


def parse_file(input: str, output: str):
    with open(input, 'r') as file:
        lines = file.readlines()
        mean_coverage = float(lines[7].split()[1])

    with open(output, 'w') as file:
        json.dump({'mean_coverage': mean_coverage}, file, indent=4)


def main():
    parser = argparse.ArgumentParser(description='Parse wgsMetrics output')
    parser.add_argument('input', help='Input file')
    parser.add_argument('output', help='Output file')
    args = parser.parse_args()

    parse_file(args.input, args.output)


if __name__ == '__main__':
    main()
