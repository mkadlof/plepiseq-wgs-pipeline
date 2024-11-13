#!/usr/bin/env python3
'''
Simple script to replace low-coverage aminoacids with N in a sequence
Input is a alignment between reference (with N in place of low coverage regions) and a sequence with SNPs and SV
Elements of the latter sequence will be replaced with Ns
1st argument is the alignment
2nd arguemnt is the name of the output fasta


Uwaga skrypt podmieniony pod grype ktora ma wiele segmentow wiec alignment i maskowanie jest w python
Uwaga skrypt oryginalnie nazywal sie w repo dedykowaneym do grypy, nanopore jako get_N_INFL.py
'''

from Bio import SeqIO
import sys
import os
from typing import Dict
import numpy as np
import subprocess
import re

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

    polecenie = ('mafft --auto --quiet --inputorder tmp.fasta')
    if 'fast' in kwargs.keys():
        polecenie = ('mafft --retree 1 --quiet --inputorder tmp.fasta') # komenda z
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
    #pozbywamy sie tepowego pliku
    os.remove('tmp.fasta')

    return slownik_alignmentu


def get_fastas(input_fasta, program):
    """
    Skrypt rozdziela pliki fasta do podplikow o nazwie program_segment.fasta
    W przypadku wystapienia w nazwie segmentu symboli "." lub "/" zostana
    one zastaione "_" w NAZWIE PLIKU
    :return:
    Dwie SPAROWANE listy jedna z nazwami segmentow bez dziwnych znakow, druga lista zwiera po prostu string
     do headera fasty
    """
    lista_plikow = []
    lista_nazw = []
    record = SeqIO.parse(input_fasta, 'fasta')
    for r in record:
        segement_name = r.id
        with open(f'{program}_{segement_name.replace("/", "_").replace(".", "_")}.fasta', 'w') as file:
            file.write(f'>{r.id}_{program}\n')
            file.write(f'{str(r.seq)}\n')
        lista_plikow.append(f'{program}_{segement_name.replace("/", "_").replace(".", "_")}.fasta')
        lista_nazw.append(f"{r.id}")
    return lista_plikow, lista_nazw

def split_final_fasta(input_fasta, prefix):
    record = SeqIO.parse(input_fasta, 'fasta')
    for r in record:
        with open(f'{prefix}_{r.id.replace("/", "_").replace(".", "_")}.fasta', 'w') as file:
            file.write(f'>{r.id}\n')
            file.write(f'{str(r.seq)}\n')

if __name__ == '__main__':
    lista_pliow_ref, _ = get_fastas(input_fasta = sys.argv[1],
                                 program='masked')
    lista_pliow_program, lista_nazw_program = get_fastas(input_fasta=sys.argv[2],
                                     program = 'medaka')
    #print(lista_pliow_ref)
    #print(lista_pliow_program)
    # output_consensus_masked.fa
    with open(f'consensus.fasta', 'w') as f:
        for plik1, plik2, fasta_header in zip(lista_pliow_ref, lista_pliow_program, lista_nazw_program):
            segment = "_".join(plik2.split('.')[0].split('_')[1:]) #
            aln = align_fasta(plik1, plik2)
            ref, target = aln.values()
            sekwencja_with_n = ''
            for i in range(len(ref)):
                if ref[i] == 'n' and target[i].upper() != 'X':
                    sekwencja_with_n += 'N'
                elif ref[i] != 'n':
                    sekwencja_with_n += target[i].upper()

            sekwencja_with_n = ''.join([element for element in sekwencja_with_n if element != 'X'])
            f.write(f'>{fasta_header}_SV\n')
            f.write(f'{sekwencja_with_n}\n')
            os.remove(plik1)
            os.remove(plik2)
    split_final_fasta(input_fasta="consensus.fasta",
                      prefix="consensus")

