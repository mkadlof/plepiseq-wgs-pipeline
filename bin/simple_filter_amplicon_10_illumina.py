#!/usr/bin/env python3

"""Zmiana w stosunku do wcześniejszej wersji - odczyty pochodzące z
dwóch amplikonów połączonych, co wynika z wypadnięcia danego primera
muszą sięgać jak wcześniej od lewego primera amplikonu -2 do prawego
primera +0, jeśli wypada primer lewy 0, a w przypadku wypadnięcia
primera prawego 0 do prawego primera amplikonu +2. Ale trafiają one do
oddzielnego pliku i sa trymowane z wykorzystaniem własnego pliku BED,
w którym mamy właściwe primer-y takiego amplikonu.

Ten skrypt wymaga programu iVar 1.4.2

"""

import subprocess
import sys
import time
from collections import defaultdict

import pysam


def read_pair_generator(bam, region_string=None):
    """Generate read pairs in a BAM file or within a region string.
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


def inner_element(amplicon_start, amplicon_end, read_start, read_end):
    return read_start > amplicon_start and read_end < amplicon_end


def write_reads(initial_bam, final_bam, statystyki, primer_left_outer, primer_left_inner,
                primer_right_outer, primer_right_inner, amplikon_bed, primer_middle, half,
                trim=True, ref_name='MN908947.3'):
    """Funkcja, która szuka odczytów w pliku initial_bam, które są w
    danym amplikonie definiowanym przez zmienne primer_left i
    primer_right. Obie zmienne sa dwuelementowymi listami, które w
    układzie - indeks half open czyli <), definiują zakres amplikonu.
    final_bam to nazwa pliku, gdzie zapiszemy bed-a, z ktorego korzysta
    ivar, samtools ampliconclip, statystyki to nazwa pliku gdzie,
    zapiszemy nazwy odczytów, które trafiały do pliku final oraz
    ich początki i końce.

    Args:
        initial_bam (str): ścieżka do pliku bam (posortowanego i
            zindeksowanego), który ma zmapowane odczyty
        final_bam (str): ścieżka do pliku, w którym planujemy trzymać
            odczyty należące do danego amplikonu
        statystyki (str): ścieżka do pliku w którym trzymamy, które
            odczyty sa w zakresie danego amplikonu
        primer_left_outer (list): trzyma granice primer'u lewego w
            postaci listy, to sa granice zewnętrzne amplikonu
        primer_left_inner (list): trzyma granice primer'u lewego w
            postaci listy, to sa granice wewnątrz amplikonu
        primer_right_outer (list): j. w. dla primera prawego
        primer_right_inner (list): j. w. dla primera prawego
        primer_middle (): to będzie g
        amplikon_bed (str): ścieżka do pliku, w którym składamy plik
            bed do filtrowania ivar/samtools
        amplikon_id (int): numer amplikonu, na którym pracujemy, na
            potrzeby tworzonego słownika

    Returns:
        int: funkcja zwraca ile, w którym podano ilość odczytów
            mapujących sie na dany amplikon
    """

    # otwieramy pliki z wejściem i wyjściem
    my_bam = pysam.AlignmentFile(initial_bam, "rb")
    my_bam_out = pysam.AlignmentFile(final_bam, "wb", template=my_bam)

    amplikon_count, statystyki_content = analyze_reads(half, my_bam, my_bam_out, primer_left_outer,
                                                       primer_middle, primer_right_outer)

    with open(statystyki, 'w') as f:
        f.write(''.join(statystyki_content))

    creating_bed(amplikon_bed, final_bam, my_bam_out, primer_left_inner, primer_left_outer,
                 primer_right_inner, primer_right_outer, ref_name, trim)

    return amplikon_count


def analyze_reads(half, my_bam, my_bam_out, primer_left_outer, primer_middle, primer_right_outer):
    amplikon_count = 0
    # Do tej zmiennej będziemy zapisywać statystyki, które później
    # zapiszemy do pliku ze statystykami
    statystyki_content = []
    # lecimy po parach odczytów
    for pair1, pair2 in read_pair_generator(my_bam):
        done = False
        # wyciagamy grance mapowan odczytu
        if pair1.reference_start <= pair2.reference_start:
            left_pair = pair1
            right_pair = pair2
        else:
            left_pair = pair2
            right_pair = pair1

        for indeks_start in range(len(primer_left_outer)):
            if done:
                continue
            for indeks_end in range(len(primer_right_outer)):
                if done:
                    continue
                begin_amplikon_L = primer_left_outer[indeks_start]
                end_amplikon_R = primer_right_outer[indeks_end]

                # illumina jestesmy poki co ostrzy jesli chodzi o kryterium zakresu
                if begin_amplikon_L <= left_pair.reference_start \
                        and end_amplikon_R > right_pair.reference_end \
                        and not done and amplikon_count <= 8000:

                    left_pair_zakres = set(range(left_pair.reference_start, left_pair.reference_end))
                    right_pair_zakres = set(range(right_pair.reference_start, right_pair.reference_end))

                    if half == 'left':
                        # wypadl prawy amplikon bierzemy prawa polowe zakresu lewy - middle - prawy
                        moj_amplikon = set(range(begin_amplikon_L, primer_middle))
                    else:
                        moj_amplikon = set(range(primer_middle, end_amplikon_R))

                    common_positions_left_pair = len(left_pair_zakres.intersection(moj_amplikon)) / float(len(left_pair_zakres))
                    common_positions_right_pair = len(right_pair_zakres.intersection(moj_amplikon)) / float(len(right_pair_zakres))

                    if (common_positions_left_pair >= 0.5 and common_positions_right_pair >= 0.1) \
                            or (common_positions_right_pair >= 0.5 and common_positions_left_pair >= 0.1):
                        # co najmniej 50 % overalpu danego read z amplikonem
                        my_bam_out.write(pair1)
                        my_bam_out.write(pair2)
                        done = True
                        amplikon_count += 1
                        statystyki_content.append(f"{pair1.qname}\t{left_pair.reference_start}\t{right_pair.reference_end}\ttwo_amplicons_{common_positions_left_pair}_{common_positions_right_pair}\n")
                    else:
                        statystyki_content.append(f"{pair1.qname}\t{left_pair.reference_start}\t{right_pair.reference_end}\ttwo_amplicons_overlap_below30\n")

                elif begin_amplikon_L <= left_pair.reference_start \
                        and end_amplikon_R > right_pair.reference_end \
                        and not done \
                        and amplikon_count > 4000:
                    statystyki_content.append(f"{pair1.qname}\t{left_pair.reference_start}\t{right_pair.reference_end}\t{'two_amplicons_above_cap'}\n")
                    done = True
    return amplikon_count, statystyki_content


def creating_bed(amplikon_bed, final_bam, my_bam_out, primer_left_inner, primer_left_outer, primer_right_inner, primer_right_outer, ref_name, trim):
    # tworzymy bed-a
    with open(amplikon_bed, 'w', encoding="utf-8") as f:
        if trim:
            for indeks_start in range(len(primer_left_outer)):
                f.write(f'{ref_name}\t{primer_left_outer[indeks_start]}\t{primer_left_inner[indeks_start]}\tncov-2019_1\t1\t+\n')
            for indeks_end in range(len(primer_right_outer)):
                f.write(
                    f'{ref_name}\t{primer_right_inner[indeks_end]}\t{primer_right_outer[indeks_end]}\tncov-2019_1\t1\t-\n')
        else:
            # tworzymy dummy trim file
            f.write(f'{ref_name}\t1\t2\tncov-2019_1\t1\t+\n')
            f.write(f'{ref_name}\t30000\t31000\tncov-2019_1\t1\t-\n')
    my_bam_out.close()
    final_bam_sort = final_bam.split('.')[0]
    pysam.sort("-o", f"{final_bam_sort}_sorted.bam", final_bam)
    pysam.index(f"{final_bam_sort}_sorted.bam")
    time.sleep(2)
    polecenie = f'ivar trim -i {final_bam_sort}_sorted.bam -b {amplikon_bed} -m 50 -q 5 -e -p {final_bam_sort}_sorted_ivar'
    ala = subprocess.Popen(polecenie, shell=True, stdout=subprocess.PIPE)
    ala.wait()


def main():
    slownik_amplikonow_with_alt_outer = {}
    slownik_amplikonow_with_alt_inner = {}
    slownik_amplikonow_uzycie = {}
    slownik_amplikonow_uzycie_left = {}  # te slowniki 3maja ile razy mamy overlap na primer lewy
    slownik_amplikonow_uzycie_rigth = {}  # te slowniki 3maja ile razy mamy overlap na primer prawy

    # ten slownik jako klucz 3ma nazwe amplikonu jak wyzej np 72
    # ale jako wartosci w odroznieniu do slonika wyzej trzyma kolejny slownik
    # z dwoma kluczami "L" i "R"
    # i dopiero te podlisty 3maja nie wartosc a liste list z grnicami amplikonow w postaci list
    # np dla amplikonu 72 ktory ma dwa primery LEFT jeden w granicach 1,10 i drug i5,20
    # oraz jeden primer prawy w granicach 100-140 mamy taki uklad
    # {'72':{'LEFT': [1, 15], 'RIGHT': [140]}}
    # jak widac dajemy zewnetrzne granice primerow bo
    # ivar w przypadku braku mapowania na oba primery odrzuci ten read gdy nie dajemy flagi e
    # tak wiec chcemy tu uwzglednic ready mapujace się na jeden primer, ktore nie 'dociagnely' do konca
    bed = sys.argv[2]
    with open(bed) as f:
        for line in f:
            line = line.split()
            # pole 3 moze miec 3 albo 4 elementy
            try:
                _, numer, kierunek, alt = line[3].split('_')
            except:
                _, numer, kierunek = line[3].split('_')

            numer = int(numer)
            if numer not in slownik_amplikonow_with_alt_outer.keys():
                slownik_amplikonow_uzycie[numer] = 0
                slownik_amplikonow_with_alt_outer[numer] = {}  # Ten slownik trzyma zewnetrzne granice amplikonow
                slownik_amplikonow_with_alt_inner[numer] = {}  # Ten slownik trzyma wewnetrzne granice amplikonow

                slownik_amplikonow_with_alt_outer[numer]['LEFT'] = []
                slownik_amplikonow_with_alt_outer[numer]['RIGHT'] = []

                slownik_amplikonow_with_alt_inner[numer]['LEFT'] = []
                slownik_amplikonow_with_alt_inner[numer]['RIGHT'] = []

                slownik_amplikonow_uzycie_left[numer] = 0
                slownik_amplikonow_uzycie_rigth[numer] = 0

            if 'LEFT' == kierunek:
                slownik_amplikonow_with_alt_outer[numer]['LEFT'].append(int(line[1]))
                slownik_amplikonow_with_alt_inner[numer]['LEFT'].append(int(line[2]))
            else:
                slownik_amplikonow_with_alt_outer[numer]['RIGHT'].append(int(line[2]))
                slownik_amplikonow_with_alt_inner[numer]['RIGHT'].append(int(line[1]))

    all_read = sys.argv[1]

    samfile_all = pysam.AlignmentFile(all_read, "rb")

    pairedreads = pysam.AlignmentFile("reads_inneramplicon.bam", "wb", template=samfile_all)

    first_pass_reject = pysam.AlignmentFile("first_pass_reject.bam", "wb", template=samfile_all)

    with open('Statystyki_one_amplicon.txt', 'w') as f:
        for pair1, pair2 in read_pair_generator(samfile_all):

            reference_start = min(pair1.reference_start, pair2.reference_start)
            reference_end = max(pair2.reference_end, pair1.reference_end)

            done = False
            for klucz in slownik_amplikonow_with_alt_outer:
                if done:
                    continue
                for indeks_start in range(len(slownik_amplikonow_with_alt_outer[klucz]['LEFT'])):
                    for indeks_end in range(len(slownik_amplikonow_with_alt_outer[klucz]['RIGHT'])):

                        begin_amplikon_L = slownik_amplikonow_with_alt_outer[klucz]['LEFT'][indeks_start]
                        end_amplikon_L = slownik_amplikonow_with_alt_inner[klucz]['LEFT'][indeks_start]

                        begin_amplikon_R = slownik_amplikonow_with_alt_inner[klucz]['RIGHT'][indeks_end]
                        end_amplikon_R = slownik_amplikonow_with_alt_outer[klucz]['RIGHT'][indeks_end]

                        # illumina jestesmy poki co ostrzy jesli chodzi o kryterium zakresu

                        if begin_amplikon_L <= reference_start \
                                and end_amplikon_R > reference_end \
                                and not done \
                                and slownik_amplikonow_uzycie[klucz] <= 4000:
                            # read jest wewnatrz amplikonu z tolerancja +2 na zakres
                            # wartosci 800 z tylnej czesci ciala, ale generalnie dzialaly na samplu 8
                            # aby zlapal wszystkie delecje
                            if begin_amplikon_L <= reference_start < end_amplikon_L:
                                slownik_amplikonow_uzycie_left[klucz] += 1
                            elif end_amplikon_R > reference_end >= begin_amplikon_R:
                                slownik_amplikonow_uzycie_rigth[klucz] += 1

                            pairedreads.write(pair1)
                            pairedreads.write(pair2)
                            done = True
                            slownik_amplikonow_uzycie[klucz] += 1
                            f.write(f"{pair1.qname}\t{reference_start}\t{reference_end}\t{'inside_amplicon'}\n")
                        elif (begin_amplikon_L + 150) <= reference_start \
                                and end_amplikon_R > reference_end \
                                and not done \
                                and slownik_amplikonow_uzycie[klucz] <= 6000:

                            if end_amplikon_R > reference_end >= begin_amplikon_R:
                                slownik_amplikonow_uzycie_rigth[klucz] += 1

                                # mam nadzieje ze w ten sposob dosampluje konce amplikonów które mogłby
                                # cierpiec z powodu bardzo duzej ilosci reaow mapujacych sie na sam poczatek amplikonu
                                # granica +150 jest wstawiona mocno arbitralnie, mozna by sie zastanowic nad przesuwaniem sie wzdluz
                                # amplikonu co jakas wartosc
                            pairedreads.write(pair1)
                            pairedreads.write(pair2)
                            done = True
                            slownik_amplikonow_uzycie[klucz] += 1
                            f.write(f"{pair1.qname}\t{reference_start}\t{reference_end}\t{'inside_amplicon_pushed'}\n")
                        elif (begin_amplikon_L + 250) <= reference_start \
                                and end_amplikon_R > reference_end \
                                and not done \
                                and slownik_amplikonow_uzycie[klucz] <= 11000:

                            if end_amplikon_R > reference_end >= begin_amplikon_R:
                                slownik_amplikonow_uzycie_rigth[klucz] += 1

                            pairedreads.write(pair1)
                            pairedreads.write(pair2)
                            done = True
                            slownik_amplikonow_uzycie[klucz] += 1
                            f.write(f"{pair1.qname}\t{reference_start}\t{reference_end}\t{'inside_amplicon_pushed_more'}\n")
                        elif begin_amplikon_L <= reference_start \
                                and end_amplikon_R > reference_end \
                                and not done \
                                and slownik_amplikonow_uzycie[klucz] > 4000:
                            f.write(f"{pair1.qname}\t{reference_start}\t{reference_end}\t{'inside_amplicon_above_cap'}\n")
                            done = True

            if not done:
                # probowalem read nie trafil do wynikowego skryptu
                # dajemy mu 2-ga szanse w kolejnym pass
                first_pass_reject.write(pair1)
                first_pass_reject.write(pair2)

        first_pass_reject.close()

        print(slownik_amplikonow_uzycie)
        print('Uzycia z lewej')
        print(slownik_amplikonow_uzycie_left)
        print('Uzycia z prawej')
        print(slownik_amplikonow_uzycie_rigth)
        print('------')
        print('Second_pass')
        pysam.sort("-o", "first_pass_reject_sorted.bam", "first_pass_reject.bam")
        pysam.index("first_pass_reject_sorted.bam")

        # drugi pass wywalamy z first_pass_reject_sorted.ba rzeczy których mapuja sie na jeden z primerow, ale jest poza drugim primerem
        # albo wypadl ten primer, albo jest to jakis artefakt, jesli jedenak primer wypadl i powstal amplikon
        # dluzszy niz powinien to takiego readu nie usuwamy (granica jest w sumie 100 odczytow mapujacych sie na amplikon),
        # jesli jest ich wiecej uznajemy ze jest to smiec i nie ma powodu go zachowywac
        # teoretycznie taki odczyt moze pochodzic z overlapujacego sie primeru z drugiego pool-a,
        # no ale skoro read startuje akurat w primer to raczej mala szansa

        samfile_second = pysam.AlignmentFile('first_pass_reject_sorted.bam', "rb")
        second_pass_reject = pysam.AlignmentFile("second_pass_reject.bam", "wb", template=samfile_second)

        for pair1, pair2 in read_pair_generator(samfile_second):
            reference_start = min(pair1.reference_start, pair2.reference_start)
            reference_end = max(pair2.reference_end, pair1.reference_end)

            done = False
            for klucz in slownik_amplikonow_with_alt_outer:
                if done:
                    continue
                for indeks_start in range(len(slownik_amplikonow_with_alt_outer[klucz]['LEFT'])):
                    for indeks_end in range(len(slownik_amplikonow_with_alt_outer[klucz]['RIGHT'])):

                        begin_amplikon_L = slownik_amplikonow_with_alt_outer[klucz]['LEFT'][indeks_start]
                        end_amplikon_L = slownik_amplikonow_with_alt_inner[klucz]['LEFT'][indeks_start]

                        begin_amplikon_R = slownik_amplikonow_with_alt_inner[klucz]['RIGHT'][indeks_end]
                        end_amplikon_R = slownik_amplikonow_with_alt_outer[klucz]['RIGHT'][indeks_end]

                        # illumina jestesmy poki co ostrzy jesli chodzi o kryterium zakresu

                        # Tak zaczałem sięidealnie w amplikonie R, ale przeskoczyłem lewy amplikon powiedzmy o 10, nie jest to
                        # read ktory bedzie w innym pool ale moge chciec takieo readu uzyc nizej wiec filtruje go
                        # tylko jesli ten amplikon ma dobre pokrycia
                        if (begin_amplikon_L - 10) <= reference_start \
                                and end_amplikon_R >= reference_end \
                                and not done and slownik_amplikonow_uzycie[klucz] > 100:

                            done = True
                            f.write(f"{pair1.qname}\t{reference_start}\t{reference_end}\t{'inside_amplicon_left_push'}\n")
                        elif begin_amplikon_L <= reference_start \
                                and (end_amplikon_R + 9) >= reference_end \
                                and not done and slownik_amplikonow_uzycie[klucz] > 100:
                            done = True
                            f.write(f"{pair1.qname}\t{reference_start}\t{reference_end}\t{'inside_amplicon_right_push'}\n")
            if not done:
                second_pass_reject.write(pair1)
                second_pass_reject.write(pair2)

        second_pass_reject.close()
        pysam.sort("-o", "second_pass_reject_sorted.bam", "second_pass_reject.bam")
        pysam.index("second_pass_reject_sorted.bam")

        for klucz, wartosc in slownik_amplikonow_uzycie.items():
            if wartosc <= 100:

                initial_bam = "second_pass_reject_sorted.bam"
                if (klucz + 2) in slownik_amplikonow_uzycie_rigth.keys():
                    # wypada prawy primer

                    final_bam = f'reads_two_amplicons_for_{klucz}_left.bam'
                    statystyki = f'Statystyki_two_amplicons_{klucz}_left.txt'
                    amplikon_bed = f'to_ivar_tmp_{klucz}_left.bed'

                    print(f'W amplikon {klucz} wypadly primer prawy i jest amplikon po prawej')
                    dodano_l = write_reads(initial_bam,
                                           final_bam,
                                           statystyki,
                                           slownik_amplikonow_with_alt_outer[(klucz)]['LEFT'],
                                           slownik_amplikonow_with_alt_inner[(klucz)]['LEFT'],
                                           slownik_amplikonow_with_alt_outer[(klucz + 2)]['RIGHT'],
                                           slownik_amplikonow_with_alt_inner[(klucz + 2)]['RIGHT'],
                                           amplikon_bed,
                                           slownik_amplikonow_with_alt_outer[klucz]['RIGHT'][0],
                                           'left',
                                           trim=True)
                else:
                    dodano_l = 0

                if (klucz - 2) in slownik_amplikonow_uzycie_rigth.keys():
                    # wypada prawy primer
                    print(f'W amplikon {klucz} wypadly primer lewy i jest amplikon po lewej')

                    final_bam = f'reads_two_amplicons_for_{klucz}_right.bam'
                    statystyki = f'Statystyki_two_amplicons_{klucz}_right.txt'
                    amplikon_bed = f'to_ivar_tmp_{klucz}_right.bed'
                    dodano_r = write_reads(initial_bam, final_bam,
                                           statystyki, slownik_amplikonow_with_alt_outer[(klucz - 2)]['LEFT'],
                                           slownik_amplikonow_with_alt_inner[(klucz - 2)]['LEFT'],
                                           slownik_amplikonow_with_alt_outer[(klucz)]['RIGHT'],
                                           slownik_amplikonow_with_alt_inner[(klucz)]['RIGHT'],
                                           amplikon_bed,
                                           slownik_amplikonow_with_alt_inner[(klucz)]['LEFT'][0],
                                           'right',
                                           trim=True)
                else:
                    dodano_r = 0

                slownik_amplikonow_uzycie[klucz] += dodano_l
                slownik_amplikonow_uzycie[klucz] += dodano_r

    print(slownik_amplikonow_uzycie)

if __name__ == '__main__':
    main()