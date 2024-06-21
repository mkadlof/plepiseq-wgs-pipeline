#!/usr/bin/env python3

"""
Script for filtering and downsampling Illumina data for influenza.
Initial filtering (proper pair and removing supplementary alignments)
is handled by samtools. Here, we will focus on additional filtering.
"""

import sys
from collections import defaultdict

import numpy as np
import pysam
from Bio import SeqIO


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    https://www.biostars.org/p/306041/
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def get_empty_amplikon_coverage_in_windows(dlugosci_segmentow, segment_id, szerokość_okna):
    """
    Funkcja vma na celu policzenie sredniego pokrycia w oknach dla segmentu. Funkcja zwraca slownik gdzie kluczami
    jest zakres okna np 0_99, 100_200 itd .... Warotoscia jest 0 ze wzgledu na specyfike danych grypowych,
    :param bam_file:
    :param segment_id:
    :return:
    """
    # Tworzymy okna i slownik ktory trzyma srednie pokrycie
    initial_split = list(range(0, dlugosci_segmentow[segment_id] + 1, szerokość_okna))
    slownik_pokrycia = {}
    for i in range(len(initial_split)):
        try:
            slownik_pokrycia[f"{initial_split[i]}_{initial_split[i + 1] - 1}"] = 0
        except IndexError:
            # jestesmy na ostatnim elemencie listy
            # poprawiamy ostatnia mozliwa granice tak by nie konczyla sie na "wartosc -1" a na wartosc
            # slownik_pokrycia[f"{initial_split[i - 1]}_{initial_split[i]}"] = 0
            del (slownik_pokrycia[f"{initial_split[i - 1]}_{initial_split[i] - 1}"])
            slownik_pokrycia[f"{initial_split[i - 1]}_{dlugosci_segmentow[segment_id] + 1}"] = 0
    return slownik_pokrycia


def read_genome_boundaries(fasta):
    """
    Wprawdzie w influenza tez korzystamy z primerow ale sa one wspolne dla wszystkich
    segmentow i usuwane wczesniej przez cutadapt-a. Na tym etapie chcemy po porstu dostac info jakie koordynaty
    maja segmenty np "chr1_PA1" ma [0,2000] TO JEST PRZYKLAD NIE RZECZYWISTE WARTOSCI
    :param fasta: sciezka do pliku fasta z genmem
    :return:
    """
    dlugosci_segmentow = {}
    genome = sequences = SeqIO.parse(fasta, "fasta")
    for segment in genome:
        dlugosci_segmentow[segment.id] = len(segment.seq)
    return dlugosci_segmentow


def read_amplicon_scheme_influenza(bed):
    """
    Dla grypy kazdy segment jest pojedycznym amplikonem ale do hierarchicznego
    wybierania odczytow zaczynamy od odczytow obejmujacych caly segment, potem schodzimy nizej
    :param bed:
    :return:
    """
    slownik_amplikonow_outer = {}  # slowanik z zewnetrznymi granicami amplikonu (czyl 5' primer left i 3' right)
    slownik_amplikonow_inner = {}  # to samo ale granice wewnetrzne (3' left i 5' right)
    slownik_amplikonow_uzycie = {}  # slownik z iloscia odczytow mapujacych sie wewnatrz danego amplikonu
    slownik_amplikonow_uzycie_left = {}  # slownik z iloscia odczytow mapujacych sie na primer left danego amplikonu
    slownik_amplikonow_uzycie_right = {}  # slownik z iloscia odczytow mapujacych sie na primer right danego amplikonu

    # kazdy segment genomu grypy to jeden segment, w celu ulatwienia parsowanai
    # amplikony maja nazwe {segmentID}_{kierunek}
    # nie ma laternatywnych primerow, zachodzacych itd ... to nie SARS
    # chr1_PB2        0       12      chr1_LEFT       1       +
    # chr1_PB2        2328    2341    chr1_RIGHT      1       -
    # chr2_PB1        0       12      chr2_LEFT       1       +

    with open(bed) as f:
        for line in f:
            line = line.split()
            segment = line[0]
            kierunek = line[3].split('_')[1]
            if segment not in slownik_amplikonow_outer.keys():
                slownik_amplikonow_outer[segment] = {}  # Ten slownik trzyma zewnetrzne granice amplikonow
                slownik_amplikonow_inner[segment] = {}  # Ten slownik trzyma wewnetrzne granice amplikonow

                slownik_amplikonow_outer[segment]['LEFT'] = 0
                slownik_amplikonow_outer[segment]['RIGHT'] = 0

                slownik_amplikonow_inner[segment]['LEFT'] = 0
                slownik_amplikonow_inner[segment]['RIGHT'] = 0

                slownik_amplikonow_uzycie[segment] = 0
                slownik_amplikonow_uzycie_left[segment] = 0
                slownik_amplikonow_uzycie_right[segment] = 0

            if 'LEFT' == kierunek:
                slownik_amplikonow_outer[segment]['LEFT'] = int(line[1])
                slownik_amplikonow_inner[segment]['LEFT'] = int(line[2])
            elif 'RIGHT' == kierunek:
                slownik_amplikonow_outer[segment]['RIGHT'] = int(line[2])
                slownik_amplikonow_inner[segment]['RIGHT'] = int(line[1])
            else:
                print('Nie rozpoznany kierunek')

    return slownik_amplikonow_outer, slownik_amplikonow_inner, \
        slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right


def filter_reads(initial_bam, final_bam_forward, final_bam_reverse, statystyki, min_length=90, mapq=30):
    all_reads = pysam.AlignmentFile(initial_bam, "rb")
    forward_reads = pysam.AlignmentFile(final_bam_forward, "wb", template=all_reads)
    reverse_reads = pysam.AlignmentFile(final_bam_reverse, "wb", template=all_reads)

    for pair1, pair2 in read_pair_generator(all_reads):

        pair1.reference_end - pair1.reference_start

        done = True
        if pair1.rlen < min_length or pair2.rlen < min_length:
            done = False
            statystyki.write(
                f"{pair1.qname}\t{pair1.reference_name}\t"
                f"{pair1.reference_start}\t{pair1.reference_end}\t"
                f"{pair2.reference_start}\t{pair2.reference_end}\tToo short\n")

        elif pair1.mapq < mapq or pair2.mapq < mapq:
            done = False
            statystyki.write(
                f"{pair1.qname}\t{pair1.reference_name}\t"
                f"{pair1.reference_start}\t{pair1.reference_end}\t"
                f"{pair2.reference_start}\t{pair2.reference_end}\tPoor mapq\n")
        elif (pair1.reference_end - pair1.reference_start) < (0.6 * pair1.rlen) or (pair2.reference_end - pair2.reference_start) < (0.6 * pair2.rlen):
            # jedne z odczytow pary ma alignment z genomem referencyjnym obejmujacy mniej niz 60% dlugosci danego odczytu
            done = False
            statystyki.write(
                f"{pair1.qname}\t{pair1.reference_name}\t"
                f"{pair1.reference_start}\t{pair1.reference_end}\t"
                f"{pair2.reference_start}\t{pair2.reference_end}\tShort alignment\n")

        if pair1.is_forward and done:
            forward_reads.write(pair1)
            reverse_reads.write(pair2)
            statystyki.write(
                f"{pair1.qname}\t{pair1.reference_name}\t"
                f"{pair1.reference_start}\t{pair1.reference_end}\t"
                f"{pair2.reference_start}\t{pair2.reference_end}\tSOK\n")
        elif pair2.is_forward and done:
            forward_reads.write(pair2)
            reverse_reads.write(pair1)
            statystyki.write(
                f"{pair1.qname}\t{pair1.reference_name}\t"
                f"{pair1.reference_start}\t{pair1.reference_end}\t"
                f"{pair2.reference_start}\t{pair2.reference_end}\tSOK\n")


if __name__ == '__main__':
    all_read = sys.argv[1]  # posortowany output z minimap-a
    bed = sys.argv[2]  # plik bed z primerami
    cap = int(sys.argv[3])  # cap na ilosc odczytow mapujacych sie na konkretny amplikon
    min_mapq = int(sys.argv[4])
    min_length = int(sys.argv[5])  # cap na smieci czyli ready nie obejmujace dwoch primerow
    reference_genome = sys.argv[6]
    statystyki = open('Statystyki.txt', 'w')

    # all_reads = pysam.AlignmentFile(initial_bam, "rb", require_index=False)

    # 0 Wczytanie informacji o dlugosci segmentow
    dlugosci_segmentow = read_genome_boundaries(reference_genome)

    # 1 Slownik z ampikonami i uzyciami amplkionow
    slownik_amplikonow_outer, slownik_amplikonow_inner, \
        slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, \
        slownik_amplikonow_uzycie_right = read_amplicon_scheme_influenza(bed=bed)

    # 2 Usuwanie odczytow za krotkich lub slabo mapujacych sie

    filter_reads(initial_bam=all_read,
                 final_bam_forward='all_reads_forward_filtered.bam',
                 final_bam_reverse='all_reads_reverse_filtered.bam',
                 statystyki=statystyki,
                 min_length=min_length,
                 mapq=min_mapq)

    # profilaktycznie randomizacja rozlozenia odczytow
    pysam.collate('-o', "all_reads_forward_filtered_randomized.bam", "all_reads_forward_filtered.bam")

    window_smoothing_bam = []
    with open('Statystyki_smieci.txt', 'w') as f:
        for segment in dlugosci_segmentow:
            slownik_slave = {}
            segment_pokrycie = get_empty_amplikon_coverage_in_windows(dlugosci_segmentow=dlugosci_segmentow,
                                                                      segment_id=segment,
                                                                      szerokość_okna=100)

            indeksy_slabych_pokryc = np.where(np.array(list(segment_pokrycie.values())) < (cap))[0].tolist()
            zakresy_zlabych_okien = [[int(klucz.split("_")[0]), int(klucz.split("_")[1])] for klucz in segment_pokrycie if segment_pokrycie[klucz] < (cap)]
            if len(zakresy_zlabych_okien) == 0:
                # nie analizuj dobrych segemtow
                continue

            all_reads = pysam.AlignmentFile("all_reads_forward_filtered_randomized.bam", "rb", require_index=False)
            pass_reads = pysam.AlignmentFile(f'{segment}_ekstra.bam', "wb", template=all_reads)

            i = 0
            for read in all_reads:
                # print(read.reference_name, read.qname)
                if read.reference_name != segment:
                    continue

                # Niech wyrownanie pokrycia idzie po read ktory jest froward
                # read 2 jest tylko slavem ktory na koniec zostanie do pliku

                reference_start = read.reference_start
                reference_end = read.reference_end
                reference_name = read.reference_name
                moj_odczyt_zakres = set(range(reference_start, reference_end))
                done = False
                # patrzymy na jakie sa pokrycia w oknach ktore obejmuje read jesli mapuje sie na slabe okno
                # done zamieniamy na True, i uzpelniamy pokrycia w slowniku
                if i % 2 != 0:
                    zakresy_zlabych_okien.reverse()

                for poor_window_start, poor_window_end in zakresy_zlabych_okien:
                    okno_zakres = set(range(poor_window_start, poor_window_end))
                    common_positions = len(okno_zakres.intersection(moj_odczyt_zakres)) / float(len(okno_zakres))
                    # print(read.qnamem, read.reference_name, reference_start, reference_end, common_positions)
                    if common_positions > 0:
                        f.write(f'Odczyt {read.qname} pochodzi z segmentu {segment} mapuje sie na {reference_start} do {reference_end} mapuje sie na okno od {poor_window_start} do {poor_window_end} i obejmuje jego {common_positions} i trafia do output\n')
                        # odczyt mapuje sine na oknx`o o slabym pokryciu
                        if not done:
                            done = True
                        # uzupelniamy slownik pokrycia
                        segment_pokrycie[f'{poor_window_start}_{poor_window_end}'] += common_positions

                # updatujemy co jest slabym oknem aby niepotrzebnie nie leciec petla przy kolejnym read
                zakresy_zlabych_okien = [[int(klucz.split("_")[0]), int(klucz.split("_")[1])] for klucz in segment_pokrycie
                                         if segment_pokrycie[klucz] < (cap)]
                if done:
                    pass_reads.write(read)
                    slownik_slave[read.qname] = 0
                else:
                    f.write(
                        f'Odczyt {read.qname} pochodzi z segmentu {segment} zakres mapowan to {reference_start} {reference_end} nie mapuje sie na okno o niskim uzyiu i nie trafia do output\n')
                i += 1

            # dodajem na koniec sparowany odczyt dla wybranych odczytow
            all_reads = pysam.AlignmentFile("all_reads_reverse_filtered.bam", "rb", require_index=False)
            for read in all_reads:
                if read.qname in slownik_slave.keys():
                    pass_reads.write(read)

            pass_reads.close()
            pysam.sort('-o', f'{segment}_ekstra_sorted.bam', f'{segment}_ekstra.bam')
            pysam.index(f'{segment}_ekstra_sorted.bam')
            window_smoothing_bam.append(f'{segment}_ekstra_sorted.bam')
    # laczenie semgentow

    pysam.merge('to_clip.bam', *window_smoothing_bam)
    pysam.sort('-o', 'to_clip_sorted.bam', 'to_clip.bam')
    pysam.index('to_clip_sorted.bam')
