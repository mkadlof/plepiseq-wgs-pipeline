#!/usr/bin/env python3

"""
Simple script to replace low-coverage aminoacids with N in a sequence
Input is a alignment between reference (with N in place of low coverage regions) and a sequence with SNPs and SV
Elements of the latter sequence will be replaced with Ns

Uwaga skrypt podmieniony pod grype ktora ma wiele segmentow wiec alignment i maskowanie jest w python
"""

import os
import re
import subprocess
import sys
from typing import Dict

import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo


def create_consensus(input_aln, segment):
    alignment = AlignIO.read(input_aln, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)
    with open(f'output_{segment}_consensus.fasta', 'w') as f:
        f.write(f'>{segment}\n')
        f.write(str(summary_align.gap_consensus(threshold=0.6, ambiguous='X')))
    return True


def align_fasta(fasta1_file: str, *args, **kwargs) -> Dict:
    """
    funkcja do alignowania sekwencji.  Sekwencje albo sa w oddzielnych plikach i wtedy laczymy je w jednym
    tmp file (usuwanym po wykonaniu funkcji), albo w jednym pliku. Do alignowania uzywamy maffta w strategii
    auto
    :param fasta1_file: str, sciezka do pliku fasta zawierajacego jedna lub wiecej sekwencji
    :param args: tupple, sciezka do dowolnej liczy plikow fasta
    :return: dict, slownik ktora jako klucze zawiera identyfikatory sekwencji z podanych pliku/plikow a jako klucze
    alignowano sekwencje
    """

    # Najpierw normalizujemy sciezki  isprawdzamy czy podane plik/pliki istnieja
    fasta1_file = [os.path.abspath(fasta1_file)]
    stan1 = os.path.isfile(fasta1_file[0])
    if not stan1:
        raise Exception('Nie podano poprawnej sciezki do pliku!')
    if len(args) > 0:
        fasta1_file.extend([os.path.abspath(x) for x in args])
        stan2 = np.alltrue([os.path.isfile(plik) for plik in fasta1_file])
        if not stan2:
            raise Exception('Podane dodatkowe pliki nie istnieja')

    for plik in fasta1_file:
        with open('tmp.fasta', 'a+') as f, open(plik, 'r') as f1:
            for line in f1:
                f.write(line)

    polecenie = ('mafft --globalpair --maxiterate 1000 --quiet --inputorder tmp.fasta')
    if 'fast' in kwargs.keys():
        polecenie = ('mafft --retree 1 --quiet --inputorder tmp.fasta')  # komenda z
    alignment = subprocess.Popen(polecenie, shell=True, stdout=subprocess.PIPE)

    # generowanie slownika
    slownik_alignmentu = {}
    for line in alignment.stdout:
        line = line.decode('UTF-8').rstrip()
        if len(re.findall(">", line)) > 0:
            klucz = line[:]
            if klucz not in slownik_alignmentu.keys():
                slownik_alignmentu[klucz] = ''
        else:
            slownik_alignmentu[klucz] = slownik_alignmentu[klucz] + line
    # pozbywamy sie tepowego pliku
    os.remove('tmp.fasta')

    return slownik_alignmentu


def get_fastas(input_fasta, program):
    """
    Skrypt rozdziela dwa pliki fasta do podplikow o nazwie program_segment.fasta
    :return:
    """
    lista_plikow = []
    record = SeqIO.parse(input_fasta, 'fasta')
    for r in record:
        with open(f'{program}_{r.id}.fasta', 'w') as file:
            file.write(f'>{r.id}_{program}\n')
            file.write(f'{str(r.seq)}\n')
        lista_plikow.append(f'{program}_{r.id}.fasta')
    return lista_plikow


def split_final_fasta(input_fasta, prefix):
    record = SeqIO.parse(input_fasta, 'fasta')
    for r in record:
        with open(f'{prefix}_{r.id}.fasta', 'w') as file:
            file.write(f'>{r.id}\n')
            file.write(f'{str(r.seq)}\n')


if __name__ == '__main__':
    lista_plikow_masked = get_fastas(sys.argv[1], 'reference')
    lista_plikow_varsacan = get_fastas(sys.argv[2], 'varscan')
    lista_plikow_freebayes = get_fastas(sys.argv[3], 'freebyes')
    lista_plikow_lofreq = get_fastas(sys.argv[4], 'lofreq')

    # 1 runda tworzenie sekwencji konsensusowe
    for masked, varscan, freebayes, lofreq in zip(lista_plikow_masked, lista_plikow_varsacan, lista_plikow_freebayes, lista_plikow_lofreq):
        # workaround for the cases when there is only one field in fasta ID
        tmp = varscan.split('_')
        if len(tmp) == 3:
            segment = tmp[2].split('.')[0]
        else:
            segment = tmp[1].split('.')[0]
        # end of workaround
        slownik_segmentu = align_fasta(varscan, freebayes, lofreq)
        with open('tmp.fasta', 'w') as f:
            for klucz, wartosc in slownik_segmentu.items():
                f.write(f'{klucz}\n')
                f.write(f'{wartosc}\n')
        create_consensus('tmp.fasta', segment)
        os.remove('tmp.fasta')
        with open(f'consensus.fasta', 'a') as f:
            aln = align_fasta(masked, f'output_{segment}_consensus.fasta')
            ref, target = aln.values()
            sekwencja_with_n = ''
            for i in range(len(ref)):
                if ref[i] == 'n' and target[i].upper() != 'X':
                    sekwencja_with_n += 'N'
                elif ref[i] != 'n':
                    sekwencja_with_n += target[i].upper()

            sekwencja_with_n = ''.join([element for element in sekwencja_with_n if element != 'X'])
            f.write(f'>{segment}\n')
            f.write(f'{sekwencja_with_n}\n')
        os.remove(masked)
        os.remove(varscan)
        os.remove(freebayes)
        os.remove(lofreq)
        os.remove(f'output_{segment}_consensus.fasta')
    split_final_fasta(f'consensus.fasta', 'consensus')
