#!/usr/bin/env python3
'''
Skrypt do parsowania pliku fastqc_data.txt. skrypt zwraca nastepujace informacje
1. wynik kazdego z testow PASS/FAILED
2. czy sa pozycje w ktorych 10th percentyl byl ponizej 30, jesli tak to ktore
3. jaki byl tile z najgorsza quality + wartosc
4. Czy jest taka overexpressed sequence ktora nie jest 'No Hit'. Jrsli tak to jaka jest to sekwencja
i jaki ma procent

w pliku fastqc_data.txt kolejne testy rozdzielone sa ">>"
'''

import sys
import os
import numpy as np
from typing import Dict, List, Any

start = 0
with open(sys.argv[1]) as f:
    failed_tests = 0
    passed_tests = 0
    bad_bases = []
    bad_bases_score = []
    bad_tiles: Dict[str, List[Any]] = {}
    overexpressed_seq = ''
    overexpressed_seq_count = 0
    overexpressed_seq_source = ''
    total_sequences = 0
    sequences_length = 0
    min_quality = sys.argv[2]

    for line in f:
        line = line.split('\t')
        if ">>" in line[0] and start == 0:
            statistics_name=line[0].strip('>>')
            statistics_outcome=line[1].rsplit()[0]
            if statistics_outcome == 'pass':
                passed_tests += 1
            else:
                failed_tests += 1
            start = 1
        elif ">>" in line[0] and start == 1:
            start = 0
            print(f'{statistics_name}\t{statistics_outcome}')
            if 'Per base sequence quality' == statistics_name:
                if len(bad_bases) == 0:
                    print("\tPoor bases\tNone\n")
                    print("\tPoor bases score\tNone\n")
                else:
                    print("\tPoor bases\t" + ";".join(map(str, [a[0] for a in sorted(zip(bad_bases, bad_bases_score), key=lambda x: x[1])])))
                    print("\tPoor bases score\t" + ";".join(map(str, [a[1] for a in sorted(zip(bad_bases, bad_bases_score), key=lambda x: x[1])])))
            elif  'Per tile sequence quality' == statistics_name:
                worst_tile = sorted(bad_tiles, key=lambda x: np.mean(bad_tiles[x]))[0]
                worst_tile_score = np.min([np.mean(x) for x in bad_tiles.values()])
                if len(bad_tiles) == 0:
                    print(f'\tWorst tile:\tNone')
                    print(f'\tWorst tile score:\tNone')
                else:
                    print(f'\tWorst tile:\t{worst_tile}')
                    print(f'\tWorst tile score:\t{worst_tile_score}')
            elif 'Overrepresented sequences' == statistics_name:
                print(f'\tSequence:\t{overexpressed_seq}')
                print(f'\tCount:\t{overexpressed_seq_count}')
                print(f'\tSource:\t{overexpressed_seq_source}')
            elif 'Basic Statistics' == statistics_name:
                print(f'\tTotal Sequences:\t{total_sequences}')
                print(f'\tSequence length:\t{sequences_length}')
            continue
        elif start == 1:
            if 'Per base sequence quality' == statistics_name:
                try:
                    if float(line[5]) < min_quality:
                        bad_bases.append(line[0])
                        bad_bases_score.append(float(line[5]))
                except:
                    # to bedzie ten naglowkowy wiersz
                    pass
            elif 'Per tile sequence quality' == statistics_name:
                try:
                    float(line[2])
                    if line[0] not in bad_tiles.keys():
                        bad_tiles[line[0]] = []
                        bad_tiles[line[0]].append(float(line[2]))
                    else:
                        bad_tiles[line[0]].append(float(line[2]))
                except:
                    #ponownie naglowek
                    pass
            elif 'Overrepresented sequences' == statistics_name:
                try:
                    if float(line[1]) > overexpressed_seq_count:
                        # podmianka w stosunku do wczesniejszej wersji sekwencje "No Hit" tez dajemy
                        # choc to pewnie tylko artefakt amplikownow
                        overexpressed_seq = line[0]
                        overexpressed_seq_count = line[1]
                        overexpressed_seq_source = line[3].rstrip()
                except:
                    pass
            elif 'Basic Statistics' == statistics_name:
                if line[0] == 'Total Sequences':
                    total_sequences = line[1].rstrip()
                elif line[0] == 'Sequence length':
                    sequences_length = line[1].rstrip()

print(f'Total no. of passed tests:\t{passed_tests}')
print(f'Total no. of failed tests:\t{failed_tests}')