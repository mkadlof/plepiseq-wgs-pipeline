#!/usr/bin/env python3

"""
A simple script that returns the total score for alignments in a BAM
file. You can also provide an argument specifying which segments
(separated by commas) are considered, e.g., chr4_HA,chr1_PB2
"""

import sys

import pysam

bam_file = sys.argv[1]

try:
    segments = sys.argv[2].split(",")
except IndexError:
    segments = []
segments = sys.argv[2].split(",")
score = 0
for read in pysam.AlignmentFile(bam_file, "rb"):
    if read.reference_name in segments and len(segments) != 0:
        score += read.get_tag('AS')
    elif read.reference_name not in segments and len(segments) != 0:
        pass
    elif len(segments) == 0:
        score += read.get_tag('AS')
    else:
        raise Exception('Unknown case')
print(score)
