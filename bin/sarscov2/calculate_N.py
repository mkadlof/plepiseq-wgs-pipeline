#!/usr/bin/env python3

"""
Simple script ot report number of N's in fasta sequence
"""

from Bio import SeqIO
import sys

record = SeqIO.parse(sys.argv[1], "fasta")
out = sys.argv[2]
with open(out, 'w') as f:
    for r in record:
        i = 0  # total number

        pos = 1  # indeksowana od 1
        slownik = {}  # tutaj trzymamy ako klucz gdzie zaczyna sie ciag N i jaka ma dlugosc
        previous = 'A'
        ciag = 0

        for element in str(r.seq).upper():
            if element == 'N':
                i += 1  # zwykle liczenie N
                if previous != 'N':
                    # zaczynamy ciag, tworzymy slownik i rozpoczynamy liczenie ciagu
                    slownik[pos] = 0
                    ciag = 1
                    pos_ciag = pos
                else:
                    ciag += 1
            else:
                if previous == 'N':
                    # konczymy ciag
                    slownik[pos_ciag] = ciag
                else:
                    # nic nie robie 
                    pass
            previous = element
            pos += 1
        # jesli konczymy N-ka zrzucamy ostatni ciag
        if previous == 'N':
            slownik[pos_ciag] = ciag

        dlugosc = len(str(r.seq))
        f.write(f'Total length of a sequence {r.id} is {len(str(r.seq))}; total numer of "N" is {i} ({i / dlugosc * 100}%); longest stretch of "N" starts at {sorted(slownik.items(), key=lambda x: x[1], reverse=True)[0][0]} and its length is {sorted(slownik.items(), key=lambda x: x[1], reverse=True)[0][1]}\n')
