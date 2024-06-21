#!/usr/bin/env python3

"""
the script returns the relative score using Needleman-Wunsch
"""

import sys

import miniseq
from minineedle import needle

fasta = miniseq.FASTA(filename=sys.argv[1])
seq1, seq2 = fasta[0], fasta[1]

alignment1: needle.NeedlemanWunsch[str] = needle.NeedlemanWunsch(seq1, seq1)
alignment1.align()
alignment2: needle.NeedlemanWunsch[str] = needle.NeedlemanWunsch(seq1, seq2)
alignment2.align()

print(alignment2.get_score() / alignment1.get_score())
