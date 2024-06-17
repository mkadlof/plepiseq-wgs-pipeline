#!/usr/bin/env python3

'''
Simple script to generate bed-like file with information about low coverage regions
nuclotides with quality below user specified threshold are not included in coverage calculations
'''

import pysam
import sys
from Bio import SeqIO

bam_file = sys.argv[1]
qual_threshold = int(sys.argv[2])
cov_threshold = int(sys.argv[3])
genome_fasta=sys.argv[4]


def read_genome_boundaries(fasta):
    genome = SeqIO.parse(fasta, "fasta")
    for segment in genome:
        return segment.id, len(segment.seq)

genome_name, genome_length = read_genome_boundaries(genome_fasta)

samfile = pysam.AlignmentFile(bam_file, "rb")

coverage = samfile.count_coverage(contig=genome_name, start=0, stop=genome_length, quality_threshold=qual_threshold)

with open('quality_mask.bed' ,'w') as f:
    for i in range(len(coverage[0])):
        suma = coverage[0][i] + coverage[1][i]  + coverage[2][i]  + coverage[3][i]
        if suma < cov_threshold:
            f.write(f'{genome_name}\t{i}\t{i + 1}\n')
