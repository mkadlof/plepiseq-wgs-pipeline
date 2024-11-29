#!/usr/bin/env python3
"""

prosty parser puszczamy ivar-a 2 razy - pierwszy raz w opcji strict tak by zapisywal tylko
te ready ktore mapuja sie na priner
drugi raz z opcja -e czyli zapisuje wszystko robi to samo co opcja strict, ale dodatkowo ready
nie mapujace sie na primery tez sa zapisywane ale z poprawionym quality check quality check

nastepnie z opcji strict pobieramy nazwy wlasciwych odczytow i dodajemy do nich ich pary z puszczania z opcji -e

Czego nie ma tutaj to readow ktore choc sa wewnatrz konkretnego amplikonu to ani jeden koniec nie mapuje sie na
zaden z dwoch primerow

Uspojnianie z illumina ? To wprawdzie dziala. ale
1. uzywac dluzszych o 1 SARS1 primers + samtools -a na odczytach interamplikon
2. dla odczyt z fuzji ktorych uzywamy jesli amplikon wypadl to  customowy bed do filtorwania
3.
"""

import pysam
import sys
import numpy as np
from Bio import SeqIO
import os

def get_amplikon_coverage_in_windows(dlugosci_segmentow, bam_file, segment_id, szerokosc_okna):
    """
    Funkcja vma na celu policzenie sredniego pokrycia w oknach dla segmentu. Funkcja zwraca slownik gdzie kluczami
    jest zakres okna np 0_99, 100_200 itd ... a wartoscia srednie pokrycie w oknie
    :param bam_file:
    :param segment_id:
    :return:
    """
    # Tworzymy okna i slownik ktory trzyma srednie pokrycie
    initial_split = list(range(0, dlugosci_segmentow[segment_id] + 1, szerokosc_okna))
    slownik_pokrycia = {}
    for i in range(len(initial_split)):
        try:
            slownik_pokrycia[f"{initial_split[i]}_{initial_split[i + 1] - 1}"] = 0
        except:
            # jestesmy na ostatnim elemencie listy
            # poprawiamy ostatnia mozliwa granice tak by nie konczyla sie na "wartosc -1" a na wartosc
            slownik_pokrycia[f"{initial_split[i - 1]}_{initial_split[i]}"] = 0


    samfile = pysam.AlignmentFile(bam_file, "rb")
    pokrycia = samfile.count_coverage(contig=segment_id, start=0, stop=dlugosci_segmentow[segment_id], quality_threshold=1)
    wektor_pokrycia = list(np.zeros(dlugosci_segmentow[segment_id]))
    for i in range(len(pokrycia[0])):
            try:
                wektor_pokrycia[i] = pokrycia[0][i] + pokrycia[1][i]  + pokrycia[2][i]  + pokrycia[3][i]
            except:
                print(f'brak pokrycia dla {segment_id} dla pozycji {i}')

    for klucz in slownik_pokrycia:
        start=int(klucz.split('_')[0])
        end = int(klucz.split('_')[1])
        slownik_pokrycia[klucz] = np.mean(wektor_pokrycia[start:end])
    return slownik_pokrycia


def calculate_coverage_in_windows():
    pass


def coverage_smoothing():
    pass


def run_mode_single(bam_file: str, chr_id: str, cycles: int, output_bam: str) -> None:
    """Run downsampling in single reads mode."""
    last_covered = -1
    reads = list(bam_file.fetch(region=chr_id))
    number_of_reads = len(reads)
    for c in range(cycles):
        print(f"Cycle {c + 1}")
        for i in range(number_of_reads):
            if reads[i] is None:
                continue
            r1 = reads[i]
            if r1.reference_start > last_covered:
                reads[i] = None
                output_bam.write(r1)
                last_covered = r1.reference_end
        last_covered = -1
    output_bam.close()


def read_genome_boundaries(fasta):
    """
    Wprawdzie w influenza tez korzystamy z primerow ale sa one wspolne dla wszystkich
    segmentow i usuwane wczesniej przez cutadapt-a. Na tym etapie chcemy po porstu dostac info jakie koordynaty
    maja segmenty np "chr1_PA1" ma [0,2000] TO JEST PRZYKLAD NIE RZECZYWISTE WARTOSCI
    :param fasta: sciezka do pliku fasta z genmem
    :return:
    """
    dlugosci_segmentow = {}
    genome = SeqIO.parse(fasta, "fasta")
    for segment in genome:
        dlugosci_segmentow[segment.id] = len(segment.seq)
    return dlugosci_segmentow

def read_amplicon_scheme_influenza(bed, bed_offset = 0):
    """
    Dla grypy kazdy segment jest pojedycznym amplikonem ale do hierarchicznego
    wybierania odczytow zaczynamy od odczytow obejmujacych caly segment, potem schodzimy nizej
    W przypadku grypy nie zawsze bedziemy mieli sekwencje z primerami (ludzie czesto przycinaja
    ja do sekwencji kodujacej. W takim wypadku bed ma strukture gdzie 2 i 3 kolumna maja identyczne
    wartosci wiec aby prefereowac dlugie odczyty sztucznie dodamy wartosc wewnetrzna)

    :param bed: standardowy plik bed z primerami
    :param bed_offset: o ile rozszerzam primery (w grypie tylko do WEWNATRZ)
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
    #chr1_PB2        0       12      chr1_LEFT       1       +
    #chr1_PB2        2328    2341    chr1_RIGHT      1       -
    #chr2_PB1        0       12      chr2_LEFT       1       +


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
                if int(line[1]) != int(line[2]):
                    slownik_amplikonow_outer[segment]['LEFT'] = int(line[1])
                    slownik_amplikonow_inner[segment]['LEFT'] = int(line[2]) + bed_offset
                else:
                    slownik_amplikonow_outer[segment]['LEFT'] = int(line[1])
                    slownik_amplikonow_inner[segment]['LEFT'] = int(line[2]) + 12 + bed_offset
            elif 'RIGHT' == kierunek:
                if int(line[1]) != int(line[2]):
                    slownik_amplikonow_outer[segment]['RIGHT'] = int(line[2])
                    slownik_amplikonow_inner[segment]['RIGHT'] = int(line[1]) - bed_offset
                else:
                    slownik_amplikonow_outer[segment]['RIGHT'] = int(line[2])
                    slownik_amplikonow_inner[segment]['RIGHT'] = int(line[1]) - 13 - bed_offset
            else:
                print('Nie rozpoznany kierunek')

    return slownik_amplikonow_outer, slownik_amplikonow_inner, \
        slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right


def _write_reads_strict_inner(initial_bam, final_bam, reject_bam,  statystyki, slownik_amplikonow_outer,
                      slownik_amplikonow_inner, slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left,
                      slownik_amplikonow_uzycie_right, dlugosci_segmentow, length_segmentu, cap = 1000, mapq = 30):
    """
    Funkcja analogiczna do tej z nanopore dla SARS-a, ale nie wymagamy posortowanego bam-a, uwzgledniamy fakt ze
    odczyty pochodza z roznych segmentow o roznej dlugosci i odfiltrowujemy nie po stalej wartosci a po
    co najmniej polowie dlugosci danego amplikonu. Pozwalamy tez na tym etapie aby odczyty byly za krotkie "do wewnatrz"
    amplikonu a nie tylko przestrzeliwaly amplikon
    :param initial_bam:
    :param final_bam:
    :param reject_bam:
    :param statystyki:
    :param slownik_amplikonow_outer:
    :param slownik_amplikonow_inner:
    :param slownik_amplikonow_uzycie:
    :param slownik_amplikonow_uzycie_left:
    :param slownik_amplikonow_uzycie_right:
    :param dlugosci_segmentow:
    :param cap:
    :param mapq:
    :return:
    """
    #
    all_reads = pysam.AlignmentFile(initial_bam, "rb", require_index=False)
    pass_reads = pysam.AlignmentFile(final_bam, "wb", template=all_reads)
    reject_reads = pysam.AlignmentFile(reject_bam, "wb", template=all_reads)
    # uwaga pysam bierze koordynaty z bam-a z 1-index na 0-index


    for read in all_reads:
        done = False

        reference_start = read.reference_start
        reference_end = read.reference_end
        reference_name = read.reference_name


        begin_amplikon_L = slownik_amplikonow_outer[reference_name]['LEFT']
        end_amplikon_L = slownik_amplikonow_inner[reference_name]['LEFT']

        begin_amplikon_R = slownik_amplikonow_inner[reference_name]['RIGHT']
        end_amplikon_R = slownik_amplikonow_outer[reference_name]['RIGHT']

        # zmiana 10.10 odczyt mapujacy sie na segment nie moze byc dluzszy niz 110% tego segmentu
        if begin_amplikon_L  <= reference_start < end_amplikon_L  \
                and end_amplikon_R  > reference_end >= begin_amplikon_R\
                and not done \
                and slownik_amplikonow_uzycie[reference_name] <= cap and read.mapq >= mapq and (reference_end - reference_start) >= (dlugosci_segmentow[reference_name] * length_segmentu)\
                and read.query_length < dlugosci_segmentow[reference_name] * 1.1:
            slownik_amplikonow_uzycie_left[reference_name] += 1
            slownik_amplikonow_uzycie_right[reference_name] += 1
            slownik_amplikonow_uzycie[reference_name] += 1
            pass_reads.write(read)
            done = True
            statystyki.write(f"{read.qname}\t{reference_name}\t{reference_start}\t{reference_end}\tinside_amplicon {reference_name} both primers \n")
        elif begin_amplikon_L  <= reference_start < end_amplikon_L \
                and end_amplikon_R  > reference_end >= begin_amplikon_R\
                and not done \
                and slownik_amplikonow_uzycie[reference_name] > cap and read.mapq >= mapq and (reference_end - reference_start) >= (dlugosci_segmentow[reference_name] * length_segmentu) \
                and read.query_length < dlugosci_segmentow[reference_name] * 1.1:
            statystyki.write(f"{read.qname}\t{reference_name}\t{reference_start}\t{reference_end}\tinside_amplicon {reference_name } both primers cap \n")
            done = True
        elif read.mapq < mapq:
            done = True
            statystyki.write(
                f"{read.qname}\t{reference_name}\t{reference_start}\t{reference_end}\tpoor quality\n")
        elif (reference_end - reference_start) < 150 and (reference_end - reference_start) > (read.rlen * 0.65):
            # 150 wzieto bo taki jest default z irmy
            done = True
            statystyki.write(
                f"{read.qname}\t{reference_name}\t{reference_start}\t{reference_end}\ttoo short\n")
        elif read.query_length >= dlugosci_segmentow[reference_name] * 1.1:
            done = True
            statystyki.write(
                f"{read.qname}\t{reference_name}\t{reference_start}\t{reference_end}\ttoo long\n")
        elif (reference_end - reference_start) >= 150 and (reference_end - reference_start) < (read.rlen * 0.65):
            # unikamy tez odczytow w ktorych mimo ze alignment jest dlugi 
            # to stanowi mniej niz 65 % dlugosci odczytu
            done = True
            statystyki.write(
                f"{read.qname}\t{reference_name}\t{reference_start}\t{reference_end}\tAlignment minorrity\n")

        if not done:
        # probowalem read nie trafil do wynikowego skryptu
        # dajemy mu 2-ga szanse w kolejnym pass
            reject_reads.write(read)

    pass_reads.close()
    reject_reads.close()
    all_reads.close()

    return slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right

def filter_reads(initial_bam, final_bam, dlugosci_segmentow, statystyki, mapq = 30, min_overlap = 0.6, max_overlap = 1.2,
                 min_alignment_overlap = 0.65, ):
    """
    Funkcja do filtrowania danych. Standardowo mamy filtrowanie na  dlugosc (podawana jako procent dlugosci segmentu),
    odrzucamy tez odczyty za dlugie. Odrzucamy zle mapujace sie odczyty oraz odczyty ktore alignuja sie tylko
    swoim fragmentem (tzn odczyt ma dlugosc 1000 ale alignment obejmuj np mniej niz 60% dlugosci takiego odczytu).
    Koncowy plik jest zrandomizowany !!!!
    :param initial_bam: Plik bam ze zmapowanymi odczytami
    :param final_bam: Plik bam tu zapisujemy output
    :param dlugosci_segmentow: slownik w ktorym kluczem jest nazwa segmentu a wartoscia jego dlugosc
    :param statystyki: otwarte lacze do pliku gdzie zapiujemy przydatne informacje
    :param mapq: int, minimalna jakosc mapowania
    :param min_overlap: jaka minimalna dlugosc musi miec odczyt (wyrazany jako procent dlugosci segmentu na jaki sie mapuje)
    :param max_overlap: jaka maksymalna dlugosc moze miec odczyt (wyrazany jako procent dlugosci segmentu)
    :param min_alignment_overlap: minimalna dlugosc alignmentu odczytu (wyrazana jako procent dlugosci odczytu)
    :return:
    """

    all_reads = pysam.AlignmentFile(initial_bam, "rb")
    all_reads_pass = pysam.AlignmentFile('tmp.bam', "wb", template=all_reads)


    for read in all_reads:
        done = True
        reference_start = read.reference_start
        reference_end = read.reference_end
        reference_name = read.reference_name
        dlugosc_segmentu = dlugosci_segmentow[reference_name]


        if read.mapq < mapq:
            done = False
            statystyki.write(f"{read.qname}\tPoor mapq\n")
        elif (reference_end - reference_start) < (read.rlen * min_alignment_overlap):
            done = False
            statystyki.write(f"{read.qname}\tAligned region is too short\n")
        elif read.rlen >  dlugosc_segmentu * max_overlap:
            done = False
            statystyki.write(f"{read.qname}\tRead is too long {read.rlen}, segment has {dlugosc_segmentu}\n")
        elif read.rlen < dlugosc_segmentu * min_overlap:
            done = False
            statystyki.write(f"{read.qname}\tRead is too long {read.rlen}, segment has {dlugosc_segmentu}\n")


        if done:
            statystyki.write(f"{read.qname}\tPassed QC\n")
            all_reads_pass.write(read)

    all_reads_pass.close()
    pysam.collate('-o', final_bam, 'tmp.bam')
    os.remove('tmp.bam')

    return True

def write_reads_strict_inner(initial_bam, final_bam, reject_bam,  statystyki, slownik_amplikonow_outer,
                      slownik_amplikonow_inner, slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left,
                      slownik_amplikonow_uzycie_right, cap = 1000):
    """
    Funkcja zapisuje odczyty ktore zaczynaja sie i koncza w primerach danego segmentu az do osiagniecia cap
    :param initial_bam: Plik bam z nieposortowanymi odczytami
    :param final_bam:  Plik bam z posortowanymi odczytami spelniajacymi kryterium mapowania
    :param reject_bam:  Plik bam z nieposortowanymi odczytami ktore Nie spelniaja kryterium
    :param statystyki:  Otwarte lacze do pliku gdzie bedziemy zapisywac statystyki
    :param slownik_amplikonow_outer: Zewnetrzne granice primerow
    :param slownik_amplikonow_inner: Wewnetrzne granice primerow
    :param slownik_amplikonow_uzycie: Ilosc readow mapujacych sie na segment
    :param slownik_amplikonow_uzycie_left:  Ilosc readow mapujacy sie na lewy primer
    :param slownik_amplikonow_uzycie_right: Ilosc readow mapujacy sie na prawy primer
    :param cap: maksymalna liczba odczytow mapujacych sie na segment zapisywan do pliku final_bam
    :return:
    """
    #
    all_reads = pysam.AlignmentFile(initial_bam, "rb", require_index=False)
    pass_reads = pysam.AlignmentFile('tmp_for_final.bam', "wb", template=all_reads)
    reject_reads = pysam.AlignmentFile(reject_bam, "wb", template=all_reads)


    for read in all_reads:
        done = False

        reference_start = read.reference_start
        reference_end = read.reference_end
        reference_name = read.reference_name


        begin_amplikon_L = slownik_amplikonow_outer[reference_name]['LEFT']
        end_amplikon_L = slownik_amplikonow_inner[reference_name]['LEFT']

        begin_amplikon_R = slownik_amplikonow_inner[reference_name]['RIGHT']
        end_amplikon_R = slownik_amplikonow_outer[reference_name]['RIGHT']

        # zmiana 10.10 odczyt mapujacy sie na segment nie moze byc dluzszy niz 110% tego segmentu
        if begin_amplikon_L  <= reference_start < end_amplikon_L  \
                and end_amplikon_R  > reference_end >= begin_amplikon_R\
                and slownik_amplikonow_uzycie[reference_name] <= cap:

            slownik_amplikonow_uzycie_left[reference_name] += 1
            slownik_amplikonow_uzycie_right[reference_name] += 1
            slownik_amplikonow_uzycie[reference_name] += 1
            pass_reads.write(read)
            done = True
            statystyki.write(f"{read.qname}\t{reference_name}\t{reference_start}\t{reference_end}\t"
                             f"inside segment {reference_name} both primers \n")
        elif begin_amplikon_L  <= reference_start < end_amplikon_L  \
                and end_amplikon_R  > reference_end >= begin_amplikon_R\
                and slownik_amplikonow_uzycie[reference_name] > cap:

            statystyki.write(f"{read.qname}\t{reference_name}\t{reference_start}\t{reference_end}\t"
                             f"inside segment {reference_name} both primers above cap \n")
            done = True
        else:
            pass

        if not done:
        # probowalem read nie trafil do wynikowego skryptu
        # dajemy mu 2-ga szanse w kolejnym pass
            reject_reads.write(read)


    pass_reads.close()
    reject_reads.close()
    all_reads.close()
    pysam.sort("-o", final_bam, 'tmp_for_final.bam')
    pysam.index(final_bam)
    os.remove('tmp_for_final.bam')

    return slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right



def get_primer_usage(initial_bam, slownik_amplikonow_outer, slownik_amplikonow_inner):
    """
    Funkcja do wyciagania informacji o ilosci odczytow mapujacych sie na dany region  wgenomie
    odpowiadajacy lokalizacji primerow
    @param initial_bam: Sciezka do pliku bam do anlizy
    @type initial_bam: basestring
    @param slownik_amplikonow_outer: slownik z informacja o granicach primerow (wartosci blizsze 3' i 5' koncom)
    @type slownik_amplikonow_outer: dict
    @param slownik_amplikonow_inner: slownik z informacja o granicach primerow (wartosc idalsze od 3' i 5')
    @type slownik_amplikonow_inner: dict
    @return: slownik z informacja jakie jest uzycie danego primeru
    @rtype: dict
    """

    uzycie_primerow = {}
    all_reads = pysam.AlignmentFile(initial_bam, "rb")
    for contig in slownik_amplikonow_outer.keys():
        left_priner_start, left_primer_end =   (slownik_amplikonow_outer[contig]['LEFT'],
                                                slownik_amplikonow_inner[contig]['LEFT'])
        right_priner_start, right_primer_end = (slownik_amplikonow_inner[contig]['RIGHT'],
                                                slownik_amplikonow_outer[contig]['RIGHT'])
        reads_left_primer = all_reads.fetch(contig=contig, start=left_priner_start, stop=left_primer_end)
        reads_right_primer = all_reads.fetch(contig=contig, start=right_priner_start, stop=right_primer_end)
        primer_usage = 0
        for x in reads_left_primer:
            primer_usage += 1
        for x in reads_right_primer:
            primer_usage += 1

        uzycie_primerow[contig] = primer_usage

    return uzycie_primerow

if __name__ == '__main__':
    all_read = sys.argv[1]  # posortowany output z minimap-a
    bed = sys.argv[2]  # plik bed z primerami
    bed_offset = int(sys.argv[3]) # rozszerz bed-a o tyle
    cap = int(sys.argv[4])  # cap na ilosc odczytow mapujacych sie na konkretny amplikon
    mapq = int(sys.argv[5])
    length_fraction_min = float(sys.argv[6]) # odczyty ktorych dlugosc alignmentu jest mniejsza niz ten procent dlugosci segmentu sa
    # odrzucane
    length_fraction_max = float(sys.argv[7])  #odczyty ktorych dlugosc alignmentu jest wieksza niz ten procent dlugosci segmentu sa
    alignment_fraction_min = float(sys.argv[8])
    reference_genome = sys.argv[9]
    window_size = int(sys.argv[10]) # wielkosc okna do wygladzania pokrycia



    statystyki = open('Statystyki.txt', 'w')

    #0 Wczytanie informacji o dlugosci segmentow
    dlugosci_segmentow = read_genome_boundaries(reference_genome)



    #1 standardowe filtrowanie
    filter_reads(initial_bam = all_read,
                 final_bam = 'all_reads_clean.bam',
                 dlugosci_segmentow = dlugosci_segmentow,
                 statystyki=statystyki,
                 min_overlap = length_fraction_min,
                 max_overlap = length_fraction_max,
                 min_alignment_overlap = alignment_fraction_min,
                 mapq = mapq)

    #1 Slownik z granicami amplikonow, ich uzyciem
    # bed offset ma tylko na celu "wzmocnienie preferencji
    # dlugich odczytow cyli takich ktore choc nie zaczynaja sie w primer
    # to sa w jego granicach
    # jesli plik bed jest dla sekwencji bez primerow ( czyli ma strukture
    # ch1_PB2 0 0 ...
    # chr1_PB2 2155 2155
    # to na sile zamienimy go na 0 12 i 2142 2155
    slownik_amplikonow_outer, slownik_amplikonow_inner, \
        slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, \
        slownik_amplikonow_uzycie_right = read_amplicon_scheme_influenza(bed=bed,
                                                                 bed_offset = bed_offset)

    slownik_uzycie_primerow = get_primer_usage(initial_bam=all_read,
                                               slownik_amplikonow_inner=slownik_amplikonow_inner,
                                               slownik_amplikonow_outer=slownik_amplikonow_outer)


    #2 W Pierwszym kroku zapisujemy odczyty ktore zaczynaja sie lub koncza w primerach (koncach segmentow)
    slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, \
        slownik_amplikonow_uzycie_right = write_reads_strict_inner(initial_bam = 'all_reads_clean.bam',
                                                                   final_bam = 'reads_inner_strict.bam',
                                                                   reject_bam = 'reject_inner_strict.bam',
                                                                   statystyki = statystyki,
                                                                   slownik_amplikonow_outer=slownik_amplikonow_outer,
                                                                   slownik_amplikonow_inner=slownik_amplikonow_inner,
                                                                   slownik_amplikonow_uzycie=slownik_amplikonow_uzycie,
                                                                   slownik_amplikonow_uzycie_left=slownik_amplikonow_uzycie_left,
                                                                   slownik_amplikonow_uzycie_right=slownik_amplikonow_uzycie_right,
                                                                   cap = cap)

    os.remove('all_reads_clean.bam')


    window_smoothing_bam = []
    with open('Statystyki_smieci.txt', 'w') as f:
        for segment in dlugosci_segmentow:
            segment_pokrycie = get_amplikon_coverage_in_windows(dlugosci_segmentow = dlugosci_segmentow,
                                                                bam_file = 'reads_inner_strict.bam',
                                                                segment_id = segment,
                                                                szerokosc_okna= window_size)
            # Czy jest takie opkno ktorego pokrycie jest mniejsze niż powiedzmy cap/4 ?
            #print(segment)
            # print(segment_pokrycie)

            indeksy_slabych_pokryc = np.where(np.array(list(segment_pokrycie.values())) < (cap))[0].tolist()

            #print(indeksy_slabych_pokryc)


            zakresy_zlabych_okien = [[int(klucz.split("_")[0]), int(klucz.split("_")[1])] for klucz in segment_pokrycie if segment_pokrycie[klucz] < (cap) ]
            #print(zakresy_zlabych_okien)
            if len(zakresy_zlabych_okien) == 0:
                # nie analizuj dobrych segemtow
                continue

            # print(f'Dla segmentu {segment} nastepujace okna maja niskie uzycie: {" ".join(map(str, indeksy_slabych_pokryc))}\n')
            # print('Odpowiada to zakresom\n')
            # print(zakresy_zlabych_okien)

            all_reads = pysam.AlignmentFile('reject_inner_strict.bam', "rb", require_index=False)
            pass_reads = pysam.AlignmentFile(f'{segment}_ekstra.bam', "wb", template=all_reads)

            # Lecimy po readach jezeli odczyt mapuje sie na słabe okno podbikamy tam uzycie o wartosc odpowiadajaca
            # overalpowi miedzy tym oknem a mapowaniem odczytu
            # uzupelniamy rowniez pokrycia w innych oknach jesli odczyt tam siega
            # i = 0

            for read in all_reads:
                #print(read.reference_name, read.qname)
                if read.reference_name != segment:
                    #print('Omijam')
                    f.write(f'Omijam odczyt {read.qname} bo nie pochodzi z segmentu {segment}\n' )
                    continue

                reference_start = read.reference_start
                reference_end = read.reference_end
                reference_name = read.reference_name
                moj_odczyt_zakres = set(range(reference_start, reference_end))
                done = False
                # patrzymy na jakie sa pokrycia w oknach ktore obejmuje read jesli mapuje sie na slabe okno
                # done zamieniamy na True, i uzpelniamy pokrycia w slowniku
                #if i % 2 != 0:
                #    zakresy_zlabych_okien.reverse()

                for poor_window_start, poor_window_end in zakresy_zlabych_okien:
                    okno_zakres=set(range(poor_window_start,poor_window_end))
                    common_positions = len(okno_zakres.intersection(moj_odczyt_zakres)) / float(len(okno_zakres))
                    #print(read.qnamem, read.reference_name, reference_start, reference_end, common_positions)
                    if common_positions > 0:
                        f.write(f'Odczyt {read.qname} pochodzi z segmentu {segment} mapuje sie na {reference_start} do {reference_end} mapuje sie na okno od {poor_window_start} do {poor_window_end} i obejmuje jego {common_positions} i trafia do output\n')
                        # odczyt mapuje sine na oknx`o o slabym pokryciu
                        if not done:
                            done = True
                        # uzupelniamy slownik pokrycia
                        segment_pokrycie[f'{poor_window_start}_{poor_window_end}'] +=  common_positions

                # updatujemy co jest slabym oknem aby niepotrzebnie nie leciec petla przy kolejnym read
                zakresy_zlabych_okien = [[int(klucz.split("_")[0]), int(klucz.split("_")[1])] for klucz in segment_pokrycie
                                         if segment_pokrycie[klucz] < (cap)]
                if done:
                    pass_reads.write(read)
                else:
                    f.write(
                        f'Odczyt {read.qname} pochodzi z segmentu {segment} zakres mapowan '
                        f'to {reference_start} {reference_end} nie mapuje sie na okno o niskim uzyiu '
                        f'i nie trafia do output\n')
                #i += 1
            pass_reads.close()
            pysam.sort('-o', f'{segment}_ekstra_sorted.bam', f'{segment}_ekstra.bam')
            pysam.index(f'{segment}_ekstra_sorted.bam')
            window_smoothing_bam.append(f'{segment}_ekstra_sorted.bam')






    # laczenie semgentow


    if len(window_smoothing_bam) > 0:
        #polacz pliki niestety pysam.merge nie zadziala bo nie wiem jak go zmusic do poprawne
        pysam.merge('window_based_smoothing.bam', *window_smoothing_bam)
        pysam.sort('-o', 'window_based_smoothing_sorted.bam', 'window_based_smoothing.bam')
        pysam.index('window_based_smoothing_sorted.bam')

        pysam.merge('to_clip.bam', 'reads_inner_strict.bam',  'window_based_smoothing_sorted.bam')
        pysam.sort('-o', 'to_clip_sorted.bam', 'to_clip.bam')
        #pysam.index('to_clip_sorted.bam')

    else:
        # wysatrcz podmienic nazwe pliku
        pysam.sort('-o', 'to_clip_sorted.bam', 'reads_inner_strict.bam')
        #pysam.index('to_clip_sorted.bam')

    with open('Primer_usage.txt', "w") as f:
        f.write('Segment\tPrimer_number\tPrimare_usage\n')
        for klucz, wartosc in slownik_uzycie_primerow.items():
            f.write(f'{klucz}\t1\t{wartosc}\n')
