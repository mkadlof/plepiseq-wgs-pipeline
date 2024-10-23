#!/usr/bin/env python3

"""
This is a simple script for generating coverage based on a quality
threshold.
"""

import sys

import pysam
from Bio import SeqIO


def read_genome_boundaries(fasta):
    """
    Although in influenza we also use primers, they are common for all
    segments and are removed earlier by cutadapt. At this stage, we
    simply want to get information about the coordinates of the
    segments, e.g., "chr1_PA1" has [0,2000]
    THIS IS AN EXAMPLE, NOT ACTUAL VALUES
    :param fasta: path to the fasta file with the genome
    :return:
    """
    segment_lengths = {}
    genome = SeqIO.parse(fasta, "fasta")
    for segment in genome:
        segment_lengths[segment.id] = len(segment.seq)
    return segment_lengths


def main():
    plik_bam = sys.argv[1]
    qual = int(sys.argv[2])
    threshold = int(sys.argv[3])
    reference_genome = sys.argv[4]

    segment_lengths = read_genome_boundaries(reference_genome)

    sam_file = pysam.AlignmentFile(plik_bam, "rb")

    for segment in segment_lengths:
        coverage = sam_file.count_coverage(contig=segment, start=0,
                                           stop=segment_lengths[segment],
                                           quality_threshold=qual)
        # coverage is a 4-element tuple with values at a given position
        # for A, C, T, and G

        with open(f'{segment.replace("/", "_")}_mask.bed', 'w') as f:
            for i in range(len(coverage[0])):
                suma = coverage[0][i] + coverage[1][i] + coverage[2][i] + coverage[3][i]
                if suma < threshold:
                    f.write(f'{segment}\t{i}\t{i + 1}\n')


if __name__ == '__main__':
    main()
