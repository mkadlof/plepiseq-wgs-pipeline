#!/usr/bin/env python3

"""
Simple script that introduces "Ns" in regions of a consensus sequnece that have low coverage
Input is an alignment between reference sequence (with N introduced in place of low coverage regions) and a consensus sequence (with introduced SNPs and SVs)
This script requires only a single argument - an alignment
"""

import sys
from Bio import SeqIO

record = list(SeqIO.parse(sys.argv[1], "fasta"))
name_of_sequence = sys.argv[1].split('.')[0].split('_')[1]

sequence_with_n = ''
i = 0  # index

for element in record[0].seq:
    if element == 'n' and record[1].seq[i] != '-':
        sequence_with_n += 'N'
    elif element != 'n':
        sequence_with_n += record[1].seq[i].upper()

    i += 1

# r remove '-' in fixed sequence
sequence_with_n = ''.join([element for element in sequence_with_n if element != '-'])
with open(f'output_{name_of_sequence}_masked.fa', 'w') as f:
    f.write(f'>{name_of_sequence}\n')
    f.write(sequence_with_n + '\n')
