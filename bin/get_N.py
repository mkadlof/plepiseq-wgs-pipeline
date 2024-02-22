#!/usr/bin/env python3

"""
Simple script to replace low-coverage aminoacids with N in a sequence
Input is a alignment between reference (with N in place of low coverage regions) and a sequence with SNPs and SV
Elements of the latter sequence will be replaced with Ns
1st argument is the alignment
2nd arguemnt is the name of the output fasta
"""

import sys

from Bio import SeqIO

record = list(SeqIO.parse(sys.argv[1], "fasta"))
name_of_variant = sys.argv[1].split('.')[0].split('_')[1]

sekwencja_with_n = ''
i = 0  # index w alignment

for element in record[0].seq:
    if element == 'n' and record[1].seq[i] != '-':
        sekwencja_with_n += 'N'
    elif element != 'n':
        sekwencja_with_n += record[1].seq[i].upper()

    i += 1

# r remove '-' in fixed sequence
sekwencja_with_n = ''.join([element for element in sekwencja_with_n if element != '-'])
with open(f'output_{name_of_variant}_masked.fa', 'w') as f:
    f.write(f'>{name_of_variant}\n')
    f.write(sekwencja_with_n + '\n')
