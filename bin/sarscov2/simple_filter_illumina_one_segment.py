#!/usr/bin/env python3

import os
import subprocess
import sys
import time
from collections import defaultdict
from typing import Dict, Any

import numpy as np
import pysam


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


def get_empty_amplikon_coverage_in_windows(slownik_amplikonow, szerokosc_okna=100):
    """
    Funkcja tworzy pusty slownik z pokryciami w oknach.
    1. Definicja zakresu amplikonu brana jest ze slownika outer, tworzonego w trakcie  wczytywania schematu amplikonow
    2. Wynikiem jest slownik slownikow. W pierwszej warstwie kluczem jest numer amlikonu i brany jest ze zmiennej
    "slownik amplikonow". Wartoscia pierwszej warstwy jest slownik w ktorej kluczem jest zakres a wartoscia int.
    Np. {1:{'0_99':0; '100_199':0}; 2:{'80_180':0; '181_280':0}}
    :param slownik_amplikonow: slownik z z zewnetrznymi granicami amplikonow, zwracany przez funkcje
    read_amplicon_scheme_RSV
    :param szerokosc_okna: szerokosc okna na jaki dzielony jest dany amplikon
    :return: opisano wyzej
    """

    slownik_pokrycia_local: Dict[Any, Dict[Any, Any]] = {}

    for klucz in slownik_amplikonow.keys():
        if klucz not in slownik_pokrycia_local:
            slownik_pokrycia_local[klucz] = {}
        # Definicja zakresu amplikonow
        start = min(slownik_amplikonow[klucz]['LEFT'])
        end = max(slownik_amplikonow[klucz]['RIGHT'])

        initial_split = list(range(start, end, szerokosc_okna))
        for i in range(len(initial_split)):
            try:
                slownik_pokrycia_local[klucz][f"{initial_split[i]}_{initial_split[i + 1] - 1}"] = 0
            except IndexError:
                # jestesmy na ostatnim elemencie listy
                # nie ma wiec wartosci o indekscie i+1
                # usuwamy poprzeni wpis i zamieniamy go tak by konczyl sie na wartosci end
                del (slownik_pokrycia_local[klucz][f"{initial_split[i - 1]}_{initial_split[i] - 1}"])
                slownik_pokrycia_local[klucz][f"{initial_split[i - 1]}_{end}"] = 0

    return slownik_pokrycia_local


def read_amplicon_scheme_RSV(bed, bed_offset=0):
    """
    Prosty skrypt do wczytywania informacji z pliku bed do slownikow.
    :param bed: Plik w formacie bed ze schematem amplikonow. Opis stosowanego formatu powinien byc umieszczony
    w dokumentacji
    :param bed_offset: int, wartosc o jaka rozszerzamy amplikon w kierunku 3' i 5' w stosunku do pliku bed
    :return: lista 4 slownikow i jeden string. slownik_amplikonow_with_alt_outer jest slownikiem slownikow. Pierwszy
    klucz prowadzi do numeru amplikonu, drugi slownik ma zawsze ten sam uklad. Dwa mozliwe klucze ('LEFT' i 'RIGHT'), a
    jako wartosci lista z zewnetrznymi granicami amplikonow. Lista jest najczesciej 1 elementowa.
    slownik_amplikonow_with_alt_outer ma identyczny uklad, ale zawiera informacje o wewnetrznych granicach amplikonow.
    slownik_amplikonow_uzycie_left i slownik_amplikonow_uzycie_right to slowniki ktore jako klucze maja numer amplikonu
    a wartosciami jest zawsze 0. 
    line[0] - nazwa genomu referencyjnego z ktorego korzystamy. Potrzebne d otworzenia tymczasowych bed-ow do ivar-a.
    """
    slownik_amplikonow_with_alt_outer: Dict[int, Dict[Any, Any]] = {}  # slowanik z zewnetrznymi granicami amplikonu
    # (czyl 5' primer left i 3' right)
    slownik_amplikonow_with_alt_inner: Dict[int, Dict[Any, Any]] = {}  # to samo ale granice wewnetrzne
    # (3' left i 5' right)
    slownik_amplikonow_uzycie_left: Dict[int, int] = {}  # slownik z iloscia odczytow mapujacych sie na primer
    # left danego amplikonu
    slownik_amplikonow_uzycie_right: Dict[int, int] = {}  # slownik z iloscia odczytow mapujacych sie na primer
    # right danego amplikonu

    with open(bed) as f:
        for line in f:
            line = line.split()
            # pole 3 moze miec 3 albo 4 elementy po split "_" bo moze wystepowac fragment "_alt"
            try:
                _, numer, kierunek, alt = line[3].split('_')
            except:
                _, numer, kierunek = line[3].split('_')
            numer = int(numer)

            if numer not in slownik_amplikonow_with_alt_outer.keys():
                slownik_amplikonow_with_alt_outer[numer] = {}  # Ten slownik trzyma zewnetrzne granice amplikonow
                slownik_amplikonow_with_alt_inner[numer] = {}  # Ten slownik trzyma wewnetrzne granice amplikonow

                slownik_amplikonow_with_alt_outer[numer]['LEFT'] = []
                slownik_amplikonow_with_alt_outer[numer]['RIGHT'] = []

                slownik_amplikonow_with_alt_inner[numer]['LEFT'] = []
                slownik_amplikonow_with_alt_inner[numer]['RIGHT'] = []

                slownik_amplikonow_uzycie_left[numer] = 0
                slownik_amplikonow_uzycie_right[numer] = 0

            if kierunek == 'LEFT':
                slownik_amplikonow_with_alt_outer[numer]['LEFT'].append(int(line[1]) - bed_offset)
                slownik_amplikonow_with_alt_inner[numer]['LEFT'].append(int(line[2]) + bed_offset)
            elif kierunek == 'RIGHT':
                slownik_amplikonow_with_alt_outer[numer]['RIGHT'].append(int(line[2]) + bed_offset)
                slownik_amplikonow_with_alt_inner[numer]['RIGHT'].append(int(line[1]) - bed_offset)
            else:
                print('Nie rozpoznany kierunek')

    return slownik_amplikonow_with_alt_outer, slownik_amplikonow_with_alt_inner, \
        slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right, line[0]


def filter_reads(initial_bam, final_bam_forward, final_bam_reverse, statystyki, min_length=90, mapq=30,
                 min_overlap=0.6):
    """
    Podstawowa funkcja do filtrowania pliku bam na dlugosc, zakres mapowania i jakosc odczytu.
    Ponadto rozbija plik fo filtrowaniu na zawierajacy mapowania 'forward'i 'reverse.
    :param initial_bam: str, sciezka do pliku bam ktory chcemy przefiltrowac
    :param final_bam_forward: str, sciezka do pliku bam, ktory ma zawierac odczyty mapujace sie na nic forward. Output
    :param final_bam_reverse: str, sciezka do pliku bam, ktory ma zawierac odczyty mapujace sie na nic reverese. Output
    :param statystyki: nazwa zmiennej z otwartym laczem do pliku do ktorego zapisywane sa informacje zwracane
    w trakcie filtrowania
    :param min_length: int, minimalna dlugosc. OBa odczyty z pary musza miec co najmniej taka dlugosc
    :param mapq: int, minimalna jakosc mapowania odczytu. OBa odczyty z pary musza miec ta wartosc.
    :param min_overlap: float, ulamek dlugosci odczytu jak mapuje sie na genom referencyjny. OBa odczyty z pary musza
     miec ta wartosc.
    :return: True, funkcja tworzy dwa pliki.
    """

    all_reads = pysam.AlignmentFile(initial_bam, "rb")
    forward_reads = pysam.AlignmentFile(final_bam_forward, "wb", template=all_reads)
    reverse_reads = pysam.AlignmentFile(final_bam_reverse, "wb", template=all_reads)

    for pair1, pair2 in read_pair_generator(all_reads):
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
        elif (pair1.reference_end - pair1.reference_start) < (0.6 * pair1.rlen) and \
                (pair2.reference_end - pair2.reference_start) < (0.6 * pair2.rlen):
            # oba z odczytow pary ma alignment z genomem referencyjnym obejmujacy mniej niz 60% dlugosci danego
            # odczytu
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
                f"{pair2.reference_start}\t{pair2.reference_end}\tFilter OK\n")
        elif pair2.is_forward and done:
            forward_reads.write(pair2)
            reverse_reads.write(pair1)
            statystyki.write(
                f"{pair1.qname}\t{pair1.reference_name}\t"
                f"{pair1.reference_start}\t{pair1.reference_end}\t"
                f"{pair2.reference_start}\t{pair2.reference_end}\tFilter OK\n")
    return True


def check_validity(slownik, cap, start, end):
    """
    Funkcja sprawdza czy dany zakres (definiowany przez argumenty start i end) zwikeszy pokrycie w regionach
    zdefiniowanych w slowniku ktorych aktualne pokrycie jest mniejsze niz cap. W takiej sytuacji zwracana jest true,
    jesli odczyt zwieksza pokrycie tylko w regionach z pokryciem juz powyzej cap to zwracany jest false
    :param slownik: slownik w ktorym kluczami sa zakresy okien (np '100_200') a jedyna wartoscia int z
     aktualnym pokryciem. To wycinek slownika tworzonego przez funkcje get_empty_amplikon_coverage_in_windows
    :param cap: int, maksymalna wartosc pokrycia w oknie
    :param start: int, poczatek zakresu
    :param end: int, koniec zakresu
    :return: bool, True jesli read jest przydatny i false jesli nie
    """
    for klucz in slownik.keys():

        klucz_start = int(klucz.split('_')[0])
        klucz_end = int(klucz.split('_')[1])

        if len(set(range(start, end)).intersection(set(range(klucz_start, klucz_end)))) > 0 and slownik[klucz] < cap:
            # read overalpuje sie na okno ktorego pokrycie jest ponizej cap
            # nie ma potrzeby sprawdzac innych okien, zwracamy True
            return True
    # odczyt nie mapuje sie na okno z niskim pokryciem
    return False


def update_slownik_pokrycia(slownik_pokrycia, reference_start, reference_end, mate_reference_start, mate_reference_end):
    """
    Funkcja do updatu wartosci w slowniku z pokryciem.
    :param slownik_pokrycia: slownik w ktorym kluczami sa zakresy okien (np '100_200') a jedyna wartoscia int z
     aktualnym pokryciem. To wycinek slownika tworzonego przez funkcje get_empty_amplikon_coverage_in_windows.
    :param reference_start: Poczatek mapowania pierwszego z odczytow z pary
    :param reference_end: Koniec mapowania pierwszego z odczytow z pary
    :param mate_reference_start: Poczatek mapowania drugiego odczytow z pary
    :param mate_reference_end: Koniec mapowania z pierwszego z odczytow z pary
    :return: update'owany slownik z pokryciem o strukturze identycznej jak ten podawany jako input
    """

    odczyt_zakres = set(range(reference_start, reference_end))
    odczyt_mate_zakres = set(range(mate_reference_start, mate_reference_end))

    for klucz in slownik_pokrycia.keys():
        # iterujemy po slowniku pokryciu klucze to zakresy np. '0_99', '100-199 itd ...

        klucz_start = int(klucz.split('_')[0])
        klucz_end = int(klucz.split('_')[1])
        okno_zakres = set(range(klucz_start, klucz_end))

        # jaka frakcja okna jest obejmowana przez odczyt i pare tego odczytu
        pokrycie_read = len(okno_zakres.intersection(odczyt_zakres)) / (klucz_end - klucz_start)
        pokrycie_mate = len(okno_zakres.intersection(odczyt_mate_zakres)) / (klucz_end - klucz_start)

        # update wartosci slownika
        slownik_pokrycia[klucz] += pokrycie_read
        slownik_pokrycia[klucz] += pokrycie_mate

    return slownik_pokrycia


def write_reads_two_amplicons(initial_bam, final_bam, statystyki, primer_left_outer, primer_left_inner,
                              primer_right_outer, primer_right_inner, amplikon_bed, primer_middle, half, initial_pokrycie, cap,
                              trim=True, ref_name='MN908947.3', lista_to_merge=[], uzyte={}):
    '''
    Tutaj mamy funkcje ktora szuka readow w pliku initial_bam, ktore sa w danym amplikonie definiowanym przez amienne
    primer_left i primer_right. Obie zmienne sa dwuelementowymi listami ktore w układzie - indek half open czyli <)
    definiuja zakres amplikonu. final_bam to nazwa pliku gdzie zapiszemy bed-a z ktorego korzysta ivar, samtools ampliconclip
    statystyki to nazwa pliku gdzie zapiszemy nazwy redow ktore trafialy do pliku final + ich poczatki i konce

    :param initial_bam: string, sciezka do pliku bam (posortowanego i zidenksowanego), ktory ma zmapowane odczyty
    :param final_bam: string sciezka do pliku w ktorym planujemy trzymac odczyty nalezace do danego amplikunu
    :param statystyki: string, sciezka do pliku w ktorym trzymamy ktore ofczyty sa w zakesie danego amplikonu
    :param primer_left_outer: lista, trzyma granice primer'u lewego w postaci listy, to sa granice zewnetrzne amplikonu
    :param primer_left_inner: lista, trzyma granice primer'u lewego w postaci listy, to sa granice wewnatrz amplikonu
    :param primer_right_outer: j.w dla primera prawrgo
    :param primer_right_inner: j.w dla primera prawego
    :param primer_middle: to bedzie granica primeru ktory wypada
    :param amplikon_bed: sciezka do pliku w ktorym skladamy plik bed do filtrowania ivarem/samtoolsami
    :param half: str, wartosc "left" lub right" w zaleznosci od tego czy amplikon o niskim pokryciu jest "z lewej" strony
    mega-amplikonu czy z prawej strony
    :param initial_pokrycie, slownik, wartosciami sa zakresy np '1_10' lub '100_200' a wartosciami inty z pokryciem
    :param cap: int, maksymalne pokrycie amplikonu
    :param trim: bool, czy przeprowadzac ad hoc maskowanie ivarem odczytow zidentyfikowanych dla tego mega amplikonu
    :param ref_name:, str, nazwa segmentu/genomu referencyjnego na ktory mapowane s odczyty. Potrzebne do tworzenia
    pliku bed do ivar-a
    :return: int, funkcja zwraca ile w ktorym podano ilosc readow mapujacych sie na dany amplikon
    '''

    # otwieramy pliki z inputem i outputem
    my_bam = pysam.AlignmentFile(initial_bam, "rb")
    my_bam_out = pysam.AlignmentFile(final_bam, "wb", template=my_bam)

    # przekazano slownik z pokryciami dla danego amplikonu
    amplikon_ilosc = initial_pokrycie

    empty = True
    with open(statystyki, 'w') as f:
        # lecimy po parach odczytów
        for pair1, pair2 in read_pair_generator(my_bam):
            done = False
            # wyciagamy grance mapowan odczytu
            # definujemy ktory z  pary odczytow mapuje sie wczesniej
            # potrzebne do okreslenia staru i konca mapowania insertu)

            if pair1.reference_start <= pair2.reference_start:
                left_pair = pair1
                right_pair = pair2
            else:
                left_pair = pair2
                right_pair = pair1

            # tutaj jak widacrozwazam oddzielnie kombinacje primerow alt i glownych
            for indeks_start in range(len(primer_left_outer)):
                if done:
                    continue
                for indeks_end in range(len(primer_right_outer)):
                    if done:
                        continue
                    begin_amplikon_L = primer_left_outer[indeks_start]
                    end_amplikon_L = primer_left_inner[indeks_start]

                    begin_amplikon_R = primer_right_inner[indeks_end]
                    end_amplikon_R = primer_right_outer[indeks_end]

                    # najpierw patrzymy czy nasz odczyt jest wewnatrz mega-amplikonu
                    # i nie przekraczamy docelowego pokrycia
                    if begin_amplikon_L <= left_pair.reference_start \
                            and end_amplikon_R > right_pair.reference_end \
                            and not done and np.any([x < cap for x in amplikon_ilosc.values()]):

                        # teraz patrzymy czy nasz odczyt wnosi cos do amplikonu ktory wypadl
                        # najpierw okreslamy zakres odczytow z pary
                        left_pair_zakres = set(range(left_pair.reference_start, left_pair.reference_end))
                        right_pair_zakres = set(range(right_pair.reference_start, right_pair.reference_end))

                        if half == 'left':
                            # wypadl prawy amplikon bierzemy lewa polowe zakresu lewy - middle - prawy
                            moj_amplikon = set(range(begin_amplikon_L, primer_middle))
                        else:
                            # wypadl lewy primer interesuje nasz prawa czesc mega-amplikonu
                            moj_amplikon = set(range(primer_middle, end_amplikon_R))

                        # okreslamy jaki procent interesujacego nas amplikonu
                        # zajmowana jest przez odczyty z pary
                        common_positions_left_pair = len(left_pair_zakres.intersection(moj_amplikon)) / float(len(left_pair_zakres))
                        common_positions_right_pair = len(right_pair_zakres.intersection(moj_amplikon)) / float(len(right_pair_zakres))

                        # TU JEST WAZNY FRAGMENT BO DECYDUJE JAKI ODCZYT JEST DODWANAY DO AMPLIKONOW
                        # O NISKIM POKRYCIU AKTUALNIE EMPIRYCZNIE WYMAGAMY BY DLA PARY JEDEN Z ODCZYTOW
                        # 50% SWOJEJ DLUGOCI OBEJMOWAL AMPLIKON O NISKIM POKRYCIU A DRUGI ODCZYT Z PARY
                        # OBEJMOWAL CO NAJMNIEJ 10% SWOJEJ DLUGOSCI TAKI AMPLIKON

                        if (common_positions_left_pair >= 0.5 and common_positions_right_pair >= 0.1) \
                                or (common_positions_right_pair >= 0.5 and common_positions_left_pair >= 0.1):
                            # co najmniej 50 % overalpu danego read z amplikonem
                            # zmienna empty posluzy nam nizej by nie robic niepotrzebnych krokow dla pustych bam-ow
                            empty = False
                            if pair1.qname not in uzyte.keys():
                                if check_validity(amplikon_ilosc, cap, left_pair.reference_start, right_pair.reference_end):
                                    # dodawaj tylko te odczyty ktore cos wnosza do pokrycia
                                    my_bam_out.write(pair1)
                                    my_bam_out.write(pair2)
                                    done = True
                                    amplikon_ilosc = update_slownik_pokrycia(amplikon_ilosc, left_pair.reference_start, left_pair.reference_end, right_pair.reference_start, right_pair.reference_end)
                                    f.write(f"{pair1.qname}\t{left_pair.reference_start}\t{right_pair.reference_end}\ttwo_amplicons_{common_positions_left_pair}_{common_positions_right_pair}\n")
                                    uzyte[pair1.qname] = 0
                        else:
                            f.write(f"{pair1.qname}\t{left_pair.reference_start}\t{right_pair.reference_end}\ttwo_amplicons_reject__{common_positions_left_pair}_{common_positions_right_pair}\n")

                    elif begin_amplikon_L <= left_pair.reference_start \
                            and end_amplikon_R > right_pair.reference_end \
                            and not done \
                            and np.all([x > cap for x in amplikon_ilosc.values()]):
                        f.write(f"{pair1.qname}\t{left_pair.reference_start}\t{right_pair.reference_end}\t{'two_amplicons_above_cap'}\n")
                        done = True

    my_bam_out.close()

    if empty:
        return amplikon_ilosc, lista_to_merge, uzyte

    # tworzymy tymczasowy plik bed ktory bedzie zawieral sie w granicach mega amplikonu
    with open(amplikon_bed, 'w') as f:
        if trim:
            for indeks_start in range(len(primer_left_outer)):
                f.write(f'{ref_name}\t{primer_left_outer[indeks_start]}\t{primer_left_inner[indeks_start]}\tdummy_1\t1\t+\n')
            for indeks_end in range(len(primer_right_outer)):
                f.write(
                    f'{ref_name}\t{primer_right_inner[indeks_end]}\t{primer_right_outer[indeks_end]}\tdummy_1\t1\t-\n')
        else:
            # nie bedzie robionego trimu, ale dla tworzenia tych samych plikow
            f.write(f'{ref_name}\t1\t2\tdummy_1\t1\t+\n')
            f.write(f'{ref_name}\t30000\t31000\tdummy_1\t1\t-\n')

    # Sortujemy plik bam
    final_bam_sort = final_bam.split('.')[0]
    pysam.sort("-o", f"{final_bam_sort}_sorted.bam", final_bam)
    pysam.index(f"{final_bam_sort}_sorted.bam")
    time.sleep(1)

    # wywolujemy ivar-a
    polecenie = (f'ivar trim -i {final_bam_sort}_sorted.bam -b {amplikon_bed} -m 30 -q 10 -e -p {final_bam_sort}_sorted_ivar')
    ala = subprocess.Popen(polecenie, shell=True, stdout=subprocess.PIPE)
    ala.wait()

    # sortowanie bam-a i indeksowanie bam
    pysam.sort("-o", f"{final_bam_sort}_sorted_ivar_sorted.bam", f"{final_bam_sort}_sorted_ivar.bam")
    pysam.index(f"{final_bam_sort}_sorted_ivar_sorted.bam")
    lista_to_merge.append(f"{final_bam_sort}_sorted_ivar_sorted.bam")

    # usuwanie plikow posrednich
    os.remove(amplikon_bed)

    os.remove(final_bam)
    os.remove(f"{final_bam_sort}_sorted.bam")
    os.remove(f"{final_bam_sort}_sorted.bam.bai")
    os.remove(f"{final_bam_sort}_sorted_ivar.bam")

    return amplikon_ilosc, lista_to_merge, uzyte


def get_inneramplicon_reads(reads_forward, reads_reverese, reads_pass, reads_reject, slownik_pokrycia,
                            slownik_amplikonow_outer, slownik_amplikonow_inner, uzycie_left, uzycie_right,
                            cap, statystyki):
    """
    Funkcja do wyrownania pokrycia tylko odczytami mapujacymi sie na jeden amplikon
    :param reads_forward: str, sciezka do pliku bam z odczytami forward, ich pozycje moga byc zrandomizowane,
    bam nie musi byc indeksowany
    :param reads_reverese: str, sciezka do pliku bam z odczytami reverese, ich pozycje moga byc zrandomizowane,
    bam nie musi byc indeksowany
    :param reads_pass: str, sciezka do pliku bam ktory zawiera odczyty forward i reverese, mapujace sie na jeden
    amplikon i ktory powinien miec w miare rowne pokrycie
    :param reads_reject: str, sciezka do pliku bam ktory zawiera odczyty forward i reverese, nie mapujace sie na
    jeden amplikon. Przydatny do analizy fuzji amplikonow
    :param slownik_pokrycia: slownik z pokryciem, wynik funkcji get_empty_amplikon_coverage_in_windows
    :param slownik_amplikonow_outer: slownik z granicami zewnetrznymi amplikonow, wynik funkcji read_amplicon_scheme_RSV
    :param slownik_amplikonow_inner: slownik z granicami wewnetrznymi amplikonow, wynik funkcji read_amplicon_scheme_RSV
    :param uzycie_left: slownik z uzyciem lewych primer'ow, wynik funkcji read_amplicon_scheme_RSV
    :param uzycie_right:, slownik z uzyciem prawych primer'ow, wynik funkcji read_amplicon_scheme_RSV
    :param cap: int, docelowe pokrycie
    :param statystyki: nazwa zmiennej z otwartym laczem do pliku do ktorego zapisywane sa informacje zwracane
    w trakcie dzialania funkcji
    :return:
    """

    all_reads_forward = pysam.AlignmentFile(reads_forward, 'rb', require_index=False)
    all_reads_reverse = pysam.AlignmentFile(reads_reverese, 'rb', require_index=False)

    pass_reads = pysam.AlignmentFile(f'{reads_pass}', 'wb', template=all_reads_forward)
    reject_reads = pysam.AlignmentFile(f'{reads_reject}', 'wb', template=all_reads_forward)

    # jako ze tylko pozycje ready forward sa randomizowane, a potrzebujemy readow reverese,
    # aby sie do nich dostac tworzymy slownik readow reverse z kluczem odpowiadajacym nazwie odczytu
    slownik_readow_reverse = {}
    for odczyt in all_reads_reverse:
        slownik_readow_reverse[odczyt.qname] = odczyt

    # idziemy po zrandomizowanych odczytach, jesli odczyt jest w amplikonie i jego para tez oraz wnosi cos do pokrycia
    # to taki odczyt zapisujemy i updateujemy slownik_pokrycia
    for odczyt in all_reads_forward:
        # gdzie mapuje sie odczyt i jego mate
        reference_start = odczyt.reference_start
        reference_end = odczyt.reference_end
        mate_reference_start = slownik_readow_reverse[odczyt.qname].reference_start
        mate_reference_end = slownik_readow_reverse[odczyt.qname].reference_end

        # okreslenie poczatkow i koncow insertu reprezentowanego przez dana pare
        pair_reference_start = min(reference_start, mate_reference_start)
        pair_reference_end = max(reference_end, mate_reference_end)

        done = False
        one_amplicon = False

        for amplikon in slownik_amplikonow_outer.keys():
            if done:
                continue
            # zewnetrzne granice amplikonu
            amplikon_start = min(slownik_amplikonow_outer[amplikon]['LEFT'])
            amplikon_end = max(slownik_amplikonow_outer[amplikon]['RIGHT'])

            # wewnetrzne granice primer'ow amplikonu, potrzebne do okreslenie uzycia primerow
            amplikon_start_inner = max(slownik_amplikonow_inner[amplikon]['LEFT'])
            amplikon_end_inner = min(slownik_amplikonow_inner[amplikon]['RIGHT'])

            if pair_reference_start >= amplikon_start and pair_reference_end < amplikon_end:
                one_amplicon = True
                # wiemy ze nasz odczy nalezy do tego amplikonu, teraz sprawdzamy
                # czy jest sens go dodawac tzn czy w oknach poprawi nam pokrycie
                if check_validity(slownik=slownik_pokrycia[amplikon],
                                  cap=cap,
                                  start=reference_start,
                                  end=reference_end):
                    # odczyt jest przydatny
                    done = True

                    # uzupelnienie slownika uzycia primerow
                    if amplikon_start_inner > pair_reference_start >= amplikon_start:
                        # odczyt zaczyna sie w lewym primer
                        uzycie_left[amplikon] += 1
                    if amplikon_end_inner <= pair_reference_end < amplikon_end:
                        # odczyt konczy sie w amplikonie prawym
                        uzycie_right[amplikon] += 1

                    # update pokrycia tym odczytem
                    slownik_pokrycia[amplikon] = update_slownik_pokrycia(slownik_pokrycia[amplikon], reference_start,
                                                                         reference_end, mate_reference_start,
                                                                         mate_reference_end)

                    # zapisanie info do pliku
                    statystyki.write(f'Odczyt {odczyt.qname} z zakresow {reference_start} {reference_end} '
                                     f'{mate_reference_end} {mate_reference_end} inside amplikonu {amplikon} below cap\n')
                else:
                    # odczyt jest nieprzydatny zapisujemy info do pliku
                    statystyki.write(f'Odczyt {odczyt.qname} z zakresow {reference_start} {reference_end} '
                                     f'{mate_reference_end} {mate_reference_end} inside amplikonu {amplikon} above cap\n')

        # zbadalem wszystkie amplikony dla danego odczytu zapisuje wynik do plikow
        if done:
            # zapisanie pary odczytow do pliku first pass
            pass_reads.write(odczyt)
            pass_reads.write(slownik_readow_reverse[odczyt.qname])
        if not one_amplicon and not done:
            # odczyt jest poza pojednyczym amplikonem i nie trafil co oczywiste do pliku first_pass

            statystyki.write(f'Odczyt {odczyt.qname} z zakresow {reference_start} {reference_end} '
                             f'{mate_reference_start} {mate_reference_end} not in amplicon\n')
            reject_reads.write(odczyt)
            reject_reads.write(slownik_readow_reverse[odczyt.qname])

    # zamkniecie plikow i ich sortowanie z indeksowaniem
    pass_reads.close()
    reject_reads.close()

    pass_sorted_name = f'{reads_pass.split(".")[0]}_sorted.bam'
    reject_sorted_name = f'{reads_reject.split(".")[0]}_sorted.bam'

    pysam.sort("-o", pass_sorted_name, reads_pass)
    pysam.index(pass_sorted_name)

    pysam.sort("-o", reject_sorted_name, reads_reject)
    pysam.index(reject_sorted_name)
    return slownik_pokrycia, uzycie_left, uzycie_right, pass_sorted_name, reject_sorted_name


if __name__ == '__main__':
    initial_bam = sys.argv[1]
    primer_bed = sys.argv[2]
    bed_offset = int(sys.argv[3])  # uwaga ivar nie poradzi sobie jesli damy cos innego niz 0
    min_length = int(sys.argv[4])
    min_mapq = int(sys.argv[5])
    window_threshold = int(sys.argv[6])
    cap = int(sys.argv[7])

    # Tutaj zapisujemy wszelkie informacje do debugu
    statystyki = open('Statystyki.txt', 'w')

    # tworzenie startowcyh slownikow z zakresami primerow oraz z uzyciem primerow
    # ivar nie poradzi sobie z poprawnym maskowanie przy bad_offset innym niz 0 wiec NA SZTYWNO ustawiam ta wartosc
    slownik_amplikonow_outer, slownik_amplikonow_inner, \
        slownik_amplikonow_uzycie_left_primer, \
        slownik_amplikonow_uzycie_right_primer, name_ref = read_amplicon_scheme_RSV(bed=primer_bed, bed_offset=0)

    # filtrowanie odczytow i rozdzielenie na ready forward i reverese
    filter_reads(initial_bam=initial_bam,
                 final_bam_forward='all_reads_forward_filtered.bam',
                 final_bam_reverse='all_reads_reverse_filtered.bam',
                 statystyki=statystyki,
                 min_length=min_length,
                 mapq=min_mapq)

    # Randomizacja odczytow forward
    pysam.collate('-o', "all_reads_forward_filtered_randomized.bam", "all_reads_forward_filtered.bam")

    # stworzenie pustego slownika z uzyciem amplikonow w oknach
    slownik_pokrycia = get_empty_amplikon_coverage_in_windows(slownik_amplikonow=slownik_amplikonow_outer,
                                                              szerokosc_okna=window_threshold)

    # Step 1 uzupelnianie pokrycia odczytami mapujacymi sie tylko na jeden konkretny amplikon. Odczyt musi byc
    # miedzy poczatkiem primeru Lewego (zamkniety) i koncem primeru prawego (otwarty)

    first_pass_name = 'first_pass.bam'
    first_reject_name = 'first_reject.bam'

    slownik_pokrycia, slownik_amplikonow_uzycie_left_primer, slownik_amplikonow_uzycie_right_primer, \
        first_pass_sorted_name, first_reject_sorted_name = \
        get_inneramplicon_reads(reads_forward='all_reads_forward_filtered_randomized.bam',
                                reads_reverese='all_reads_reverse_filtered.bam',
                                reads_pass=first_pass_name,
                                reads_reject=first_reject_name,
                                slownik_pokrycia=slownik_pokrycia,
                                slownik_amplikonow_outer=slownik_amplikonow_outer,
                                slownik_amplikonow_inner=slownik_amplikonow_inner,
                                uzycie_left=slownik_amplikonow_uzycie_left_primer,
                                uzycie_right=slownik_amplikonow_uzycie_right_primer,
                                cap=cap,
                                statystyki=statystyki)

    # Drugi krok, odsianie odczytow ktore mapuja sie przede wszystkim na "dobry mplikon" a ktore "przestrzeliły"
    # jego granice o max 10 bp. Jest to bardzo empiryczny threshold. Przestrzelony amplikon musi miec pokrycie powyzej
    # 100 aby taki odczyt nie byl zapisany. Jest t owiec forma filtra, aby sekwencje mapujace sie de facto
    # to ampliko o dobrym uzyciu Nie byly rozwazane w kroku trzeciem

    # input kroku
    # first_reject_sorted_bam = pysam.AlignmentFile(first_reject_sorted_name, 'rb', require_index=True)
    # Dodatkowe uzywane zmienne globalne dla tego kroku
    # 1. slownik_amplikonow_outer
    # 2. slownik_pokrycia
    # 3. statystyki

    # output kroku
    ## sekwencje ktore NIE sa przestrzelone o +/- 10 bp a obejmuja zakres
    ## wiekszy niz pojedynczy amplikon
    # second_pass_name = 'second_pass.bam'
    # second_pass_bam = pysam.AlignmentFile(second_pass_name, "wb", template=first_reject_sorted_bam)
    #
    # #parametry kroku
    # empirical_cap = 100
    # zachowane = 0
    #
    # # funkcja
    # for pair1, pair2 in read_pair_generator(first_reject_sorted_bam):
    #     pair_reference_start = min(pair1.reference_start, pair2.reference_start)
    #     pair_reference_end = max(pair1.reference_end, pair2.reference_end)
    #     done = False
    #
    #     for amplikon in slownik_amplikonow_outer.keys():
    #         if done:
    #             continue
    #         # zewnetrzne granice amplikonu
    #         amplikon_start = min(slownik_amplikonow_outer[amplikon]['LEFT'])
    #         amplikon_end = max(slownik_amplikonow_outer[amplikon]['RIGHT'])
    #         # srednie pokrycie w amplikonie
    #         amplikon_pokrycie = np.mean(list(slownik_pokrycia[amplikon].values()))
    #
    #         if (amplikon_start - 10) <= pair_reference_start < amplikon_start and pair_reference_end < amplikon_end \
    #             and amplikon_pokrycie > empirical_cap:
    #             done = True
    #             statystyki.write(f'Odczyt {pair1.qname} z zakresow {pair_reference_start} {pair_reference_end} '
    #                              f'przestrzelony start ale over cap w amplikon {amplikon}\n')
    #
    #         elif amplikon_start <= pair_reference_start and pair_reference_end < (amplikon_end + 10) \
    #             and amplikon_pokrycie > empirical_cap:
    #             done = True
    #             statystyki.write(f'Odczyt {pair1.qname} z zakresow {pair_reference_start} {pair_reference_end} '
    #                              f'przestrzelony end ale over cap w amplikon {amplikon}\n')
    #
    #     if not done:
    #         zachowane +=1
    #         second_pass_bam.write(pair1)
    #         second_pass_bam.write(pair2)
    #         statystyki.write(f'Odczyt {pair1.qname} z zakresow {pair_reference_start} {pair_reference_end} '
    #                          f'two amplikons\n')
    #
    #
    #
    # second_pass_bam.close()
    # second_pass_sorted_name = f'{second_pass_name.split(".")[0]}_sorted.bam'
    #
    # if zachowane > 0:
    #     pysam.sort("-o", second_pass_sorted_name, second_pass_name)
    #     pysam.index(second_pass_sorted_name)
    # else:
    #     # funkcja wyzej nie odrzucila zadnego readu wiec do analizy fucji amplikonow
    #     # wchodzi "goly" output pierwszego kroku
    #     pysam.sort("-o", second_pass_sorted_name, first_reject_sorted_name)
    #     pysam.index(second_pass_sorted_name)

    first_reject_sorted_bam = pysam.AlignmentFile(first_reject_sorted_name, 'rb', require_index=True)

    second_pass_sorted_name = f'second_pass_sorted.bam'

    pysam.sort("-o", second_pass_sorted_name, first_reject_sorted_name)
    pysam.index(second_pass_sorted_name)

    # Ostatni krok dla amplikonow o niskim pokryciu dodajemy odczyty ktore moga pochodzic z fuzji
    # z amplikonami +2 lub -2. Takie odczyty dalej musza znaczaco overlapowac amplikon o niskim uzyciu
    # co najmniej 60% jednego z odczytow z pary i 10 % drugiego z odczytow z pary

    # to merge trzyma nazwy powstalych plikow (kazdy fuzja amplikonow generuje oddzielny plik).
    # Potem te pliki laczymy i sortujemy
    to_merge = []
    # ten slownik trzyma nazwy pary odczytow ktore zostaly uzyte jako powstale z fuzji. Dzieki temu nie dodamy
    # danej pary odczytow wiecej niz raz
    uzyte_odczyty = {}

    empirical_cap = 100

    # iterujemy po slowniku pokrycia do analizy wejda tylko takie amplikony, ktorych co najmniej jedno
    # okno z zakresu amplikonu ma pokrycie ponizej cap
    for klucz, wartosc in slownik_pokrycia.items():

        # amplikon_pokrycie = np.mean(list(wartosc.values()))
        if np.any(np.array(list(wartosc.values())) < empirical_cap):
            # procesujemy jesli dowolne okno w amplikonie ma pokrycie mniejsze niz cap
            # readu poza mapowaniem na amplikon musi mapowac sie na okno w tym amplikonie
            # inaczej odczyt jest bez senwsu
            initial_bam = second_pass_sorted_name

            if (klucz + 2) in slownik_pokrycia.keys():
                # szukamy fuzji amplikonu z amplikonem + 2
                # ustawianie tymczasowych plikow
                final_bam = f'reads_two_amplicons_for_{klucz}_left.bam'
                statystyki_two_amplicons = f'Statystyki_two_amplicons_{klucz}_left.txt'
                amplikon_bed = f'to_ivar_tmp_{klucz}_left.bed'
                tekst = '\t'.join(map(str, list(wartosc.values())))
                statystyki.write(f'W amplikon {klucz} pokrycie w oknach wynioslo {tekst} '
                                 f'i jest amplikon po prawej\n')

                # slownik_pokycia_fuzja = {**slownik_pokrycia[klucz], **slownik_pokrycia[klucz + 2]}

                slownik_pokrycia[klucz], to_merge, uzyte_odczyty = write_reads_two_amplicons(
                    initial_bam=initial_bam,
                    final_bam=final_bam,
                    statystyki=statystyki_two_amplicons,
                    primer_left_outer=slownik_amplikonow_outer[(klucz)]['LEFT'],
                    primer_left_inner=slownik_amplikonow_inner[(klucz)]['LEFT'],
                    primer_right_outer=slownik_amplikonow_outer[(klucz + 2)]['RIGHT'],
                    primer_right_inner=slownik_amplikonow_inner[(klucz + 2)]['RIGHT'],
                    amplikon_bed=amplikon_bed,
                    primer_middle=max(slownik_amplikonow_outer[klucz]['RIGHT']),
                    half='left',
                    trim=True,
                    cap=empirical_cap,
                    initial_pokrycie=slownik_pokrycia[klucz],
                    ref_name=name_ref,
                    lista_to_merge=to_merge,
                    uzyte=uzyte_odczyty)

                tekst = '\t'.join(map(str, slownik_pokrycia[klucz]))
                statystyki.write(f'W amplikon {klucz} pokrycie wzroslo do {tekst} '
                                 f'po fuzji z amplikonem + 2\n')
            # Przeliczenie pokrycia po pierwszym kroku fuzjo
            amplikon_pokrycie = np.mean(list(wartosc.values()))

            if (klucz - 2) in slownik_pokrycia.keys() and np.any(np.array(list(slownik_pokrycia[klucz].values())) < empirical_cap):
                # wypada prawy primer i pokrycie jest ciagle ponizej cap

                statystyki.write(f'W amplikon {klucz} pokrycie srednie wynioslo {amplikon_pokrycie} po fuzji +2'
                                 f' jest dalej nizsze niz {empirical_cap} dodaje odczyty z fuzji z maplikonem -2\n')

                final_bam = f'reads_two_amplicons_for_{klucz}_right.bam'
                statystyki_two_amplicons = f'Statystyki_two_amplicons_{klucz}_right.txt'
                amplikon_bed = f'to_ivar_tmp_{klucz}_right.bed'

                # slownik_pokycia_fuzja = {**slownik_pokrycia[klucz - 2 ] , **slownik_pokrycia[klucz]}

                slownik_pokrycia[klucz], to_merge, uzyte_odczyty = write_reads_two_amplicons(
                    initial_bam=initial_bam,
                    final_bam=final_bam,
                    statystyki=statystyki_two_amplicons,
                    primer_left_outer=slownik_amplikonow_outer[(klucz - 2)]['LEFT'],
                    primer_left_inner=slownik_amplikonow_inner[(klucz - 2)]['LEFT'],
                    primer_right_outer=slownik_amplikonow_outer[(klucz)]['RIGHT'],
                    primer_right_inner=slownik_amplikonow_inner[(klucz)]['RIGHT'],
                    amplikon_bed=amplikon_bed,
                    primer_middle=min(slownik_amplikonow_inner[(klucz)]['LEFT']),
                    half='right',
                    trim=True,
                    cap=empirical_cap,
                    initial_pokrycie=slownik_pokrycia[klucz],
                    ref_name=name_ref,
                    lista_to_merge=to_merge,
                    uzyte=uzyte_odczyty)
                tekst = '\t'.join(map(str, slownik_pokrycia[klucz]))
                statystyki.write(f'W amplikon {klucz} pokrycie wzroslo do {tekst} '
                                 f'po fuzji z amplikonem - 2\n')

    # laczenie semgentow jesli byly bam-y zawierajace odczyty z fuzji amplikonow
    if (len(to_merge) > 0):
        pysam.merge('two_amplicons.bam', *to_merge)
        pysam.sort('-o', 'two_amplicons_sorted.bam', 'two_amplicons.bam')
        pysam.index('two_amplicons_sorted.bam')

    # Uzycie primerow
    with open('Primer_usage.txt', "w") as f:
        f.write('Segment\tPrimer_number\tPrimare_usage\n')
        for klucz in slownik_amplikonow_uzycie_left_primer.keys():
            f.write(f'{name_ref}\t{klucz}\t{int(slownik_amplikonow_uzycie_left_primer[klucz]) + int(slownik_amplikonow_uzycie_right_primer[klucz])}\n')

    # Legacy version for Statisics file
    statystyki.write('Uzycie primerow (kolejno numer primeru, uzycie primeru lewego, uzycie primeru prawego:\n')

    for klucz in slownik_amplikonow_uzycie_left_primer.keys():
        statystyki.write(f'{klucz}\t{slownik_amplikonow_uzycie_left_primer[klucz]}\t'
                         f'{slownik_amplikonow_uzycie_right_primer[klucz]}\n')
