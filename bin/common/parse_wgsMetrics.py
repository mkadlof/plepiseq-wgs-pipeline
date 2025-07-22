#!/usr/bin/env python3

import argparse
import json


def parse_file(input: str, output: str):
    # Extract mean coverage from picard stats file
    with open(input, 'r') as file:
        lines = file.readlines()
        mean_coverage = float(lines[7].split()[1])

    # Extract coverage_histogram_data from picard stats file
    with open(input, 'r') as file:
        extract = False
        lines = []
        for line in file:
            if line.startswith("## HISTOGRAM"):
                extract = True
                continue
            if extract:
                stripped_line = line.strip()
                if stripped_line:  # Omit empty lines
                    lines.append(stripped_line)
    with open("coverage_histogram_data.tsv", 'w') as file:
        file.write('\n'.join(lines))

    with open(output, 'w') as f:
        f.write(json.dump({'average_coverage_value': mean_coverage}, ensure_ascii=False, indent=4))


def main():
    parser = argparse.ArgumentParser(description='Parse wgsMetrics output')
    parser.add_argument('picard_stats', help='Input file (picard stats)')
    parser.add_argument('output', help='Output file')
    args = parser.parse_args()

    parse_file(args.picard_stats, args.output)


if __name__ == '__main__':
    main()
