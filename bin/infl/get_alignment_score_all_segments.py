#!/usr/bin/env python3

"""
A simple script that returns the total score for alignments in a BAM
file but unlike the original script, it returns the results for ALL
segments in the form of a row matrix
"""

import sys

import pysam

bam_file = sys.argv[1]

# We need to have separate lists for influenza A and B
segments = 'chr1_PB2,chr2_PB1,chr3_PA,chr4_HA,chr5_NP,chr6_NA,chr7_MP,chr8_NS'.split(",")

scores = {key: 0 for key in segments}
for read in pysam.AlignmentFile(bam_file, "rb"):
    if read.reference_name in segments:
        scores[read.reference_name] += read.get_tag('AS')
    else:
        print(f'Unknown segment for {read.qname}')

text = ''
for klucz in segments:
    text = text + f'{scores[klucz]}\t'

print(text)
