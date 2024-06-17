#!/usr/bin/env python3

"""
Another simple script to parse ivar output (which is send to out.log in the main pipeline). To determine
1. IF a given amplicon pairs does not work (including all alternative versions of L/R primer in a given region)
2. We have a potential dip in coverage due to 3 consecutives primer not workin, this is an expantion

-----<YES-----------NOT>----------------------<NOT--------YES>---------------------------  <- pool1
--------------------------------------NNNNNNNNNNNNNNNN--------------------------------------------
------------------------------<YYES---NOT>-----------------------------<YES--------------- <- pool2
"""

import sys

#parsujemy beda
# tworzymy slownik nazw primerow gdzie kluczem jest numer porzadkowy np 10
# wartosciami jest 2-elementowa lista w ktorej pierwszy element t liczba readow z Lewy primerem (glownym i wszystkimi pomocniczymi)
# a drugi element to samo ale dla primeru Prawego
slownik_uzycia = {}
# pierwszy element to bed ze schematem wszystkich amplikowno
with open(sys.argv[1]) as f1:
    for line in f1:
        line=line.split('\t')
        indeks = int(line[3].split('_')[1])
        if indeks not in slownik_uzycia:
            slownik_uzycia[indeks] = [0,0]

# drugi argument to plik z outputem gdzie raportowne jestuzycie ivara
with open(sys.argv[2]) as f1:
    for line in f1:
        if 'nCoV-2019_' in line:
            # aktualnie tylko primery maja taka nazwe
            nazwa, uzycie = line.split('\t')
            nazwa = nazwa
            uzycie = int(uzycie)
            indeks,kierunek = nazwa.split('_')[1:3] # ignorujemy czuj jest alternatywny
            indeks = int(indeks)
            if kierunek == 'LEFT':
                slownik_uzycia[indeks][0] +=  uzycie
            elif kierunek == 'RIGHT':
                slownik_uzycia[indeks][1] += uzycie
            else:
                raise Exception(f'Nieobslugiwany wyjatek dla {line}')

# Test czy sa takie indeksy ktore nie uzywne ponizej thresholdu

threshold = int(sys.argv[3])
cons_len = 0
poor_pair =[]
poor_multiple = []

for key in sorted(slownik_uzycia.keys()):
    if slownik_uzycia[key][0] < int(threshold/2):
        # albo kontynuujemy albo rozpoczynamy wstawke
        cons_len +=1
    else:
        if cons_len >= 3:
            # zamykamy dziure 3 i wiecej primerow
            poor_multiple.append(f'{cons_len}<-{key}_LEFT')
            #print(f'Poor primer usage of length {cons_len} ending prior to {key}_LEFT')
        cons_len = 0

    if slownik_uzycia[key][1] < int(threshold/2):
        cons_len += 1
    else:
        if cons_len >= 3:
            # zamykamy dziure 3 i wiecej primerow
            poor_multiple.append(f'{cons_len}<-{key}_RIGHT')
            #print(f'Poor primer usage of length {cons_len} ending prior to {key}_RIGHT')
        cons_len = 0
    #if slownik_uzycia[key][0] < threshold and slownik_uzycia[key][1] < threshold:
    if (slownik_uzycia[key][0] + slownik_uzycia[key][1]) < threshold:
        poor_pair.append(key)

if cons_len >= 3:
    poor_multiple.append(f'{cons_len}<-{key}_RIGHT')
print(f'Poor primers stretch:\t{",".join(poor_multiple)}')
print(f'Poor primers pair:\t{",".join(map(str, poor_pair ))}')
# Kiedy dochodzi do 3 nieuzywanych primerow




