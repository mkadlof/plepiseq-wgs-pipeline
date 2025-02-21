#!/usr/bin/python
# w Etoki analizie podlegaja tylko conti o dlugosci 300+
from Bio import SeqIO
import sys
import re
import numpy as np

def read_genome(fasta):
    """
    Funkcja do wczytania genomu po 3 iteracjach pilonem
    Zwracany slownik ma jako klucz nazwe genomu, a jako wartosc liste
    zawierajaca w kolejnosci dlugosc contiga, pokrycie i sekwencje
    :param fasta:
    :return:
    """
    slownik = {}
    genome = SeqIO.parse(fasta, "fasta")
    for segment in genome:
        slownik[segment.id] = []
        slownik[segment.id].append(len(segment.seq))
        slownik[segment.id].append(float(re.search(r'(\d+\.\d+)', segment.id)[0]))
        slownik[segment.id].append(str(segment.seq))
    return slownik


# nazywam to fastq bo etoki 3ma taka notacje
# jedyne czego nie mam to ciagu quality
fastq = read_genome(sys.argv[1])
fastq_with_rejects = read_genome(sys.argv[2])

# Uwaga bierzemy contigi o dlugosci > 10 (oryginalnie bylo 300 jak w enterobase). Zgrac to z Tomkiem
seq = sorted([s for s in fastq.values() if s[0] >= 300], key=lambda x: -x[0])
n_seq = len(seq)
n_base = sum([s[0] for s in seq])
n50, acc = 0, [0, 0]
l50, ave_depth, n_low = 0, 0, 0
if n_seq > 0:
    for l50, s in enumerate(seq):
        acc[0], acc[1] = acc[0] + s[0], acc[1] + s[0] * s[1]
        if acc[0] * 2 >= n_base:
            n50 = s[0]
            break
    l50 += 1
    ave_depth = acc[1] / acc[0]
    for s in seq:
        n_low += np.sum(np.array(list(s[2])) == 'N')
with open('Summary_statistics.txt', 'w') as f:
    f.write(f'n_contig = {n_seq}\nn_base = {n_base}\nave_depth = {ave_depth}\nn_N = {float(n_low)}\nN50 = {n50}\nL50 = {l50}\n')

# Powtarzamy statystyki ale rowniez dla odrzuconych alleli (i dodatkowo z minimalnych thresholdem na dlugosc)
# wyniki zostana zapisane do innego pliku nie zwyniki dla przefiltrowanego genomu
seq = sorted([s for s in fastq_with_rejects.values() if s[0] >= 10], key=lambda x: -x[0])
n_seq = len(seq)
n_base = sum([s[0] for s in seq])
n50, acc = 0, [0, 0]
l50, ave_depth, n_low = 0, 0, 0
if n_seq > 0:
    for l50, s in enumerate(seq):
        acc[0], acc[1] = acc[0] + s[0], acc[1] + s[0] * s[1]
        if acc[0] * 2 >= n_base:
            n50 = s[0]
            break
    l50 += 1
    ave_depth = acc[1] / acc[0]
    for s in seq:
        n_low += np.sum(np.array(list(s[2])) == 'N')


with open('Summary_statistics_with_reject.txt', 'w') as f:
    f.write(f'n_contig = {n_seq}\nn_base = {n_base}\nave_depth = {ave_depth}\nn_N = {float(n_low)}\nN50 = {n50}\nL50 = {l50}\n')
