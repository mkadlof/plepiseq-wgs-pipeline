#!/usr/bin/env python3

"""
Parsujemy vcf-a ktory ma zawierac ambigous nucleotides,
podmieniamy WSZYSTKIE pozycje z tego vcf-a wiec filtrowanie robimy gdzie indziej ...
"""

import sys

slownik_amb = {'R': ['A', 'G'],
               'Y': ['T' , 'C'],
               'K': ['G', 'T'],
               'M': ['A', 'C'],
               'S':['G', 'C'],
               'W':['A', 'T']}

with open(sys.argv[1]) as f1, open(sys.argv[2], 'w') as f2:
    for line in f1:
        if "#" in line:
            f2.write(line)
        else:
            # linia bez # czyli nasza mutacja
            line = line.split('\t')
            ref = line[3]
            alt = line[4]
            rest = "\t".join(line[5:])
            new_alt = [k for k, v in slownik_amb.items() if sorted([ref, alt]) == sorted(v)][0]
            f2.write(f'{line[0]}\t{line[1]}\t{line[2]}\t{ref}\t{new_alt}\t{rest}')
