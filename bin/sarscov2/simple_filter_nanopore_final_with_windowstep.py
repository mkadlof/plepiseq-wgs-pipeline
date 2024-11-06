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
import os
from typing import Dict, Any


def read_amplicon_scheme(bed, bed_offset=1):
    """
    Prosta funkcja to wczytywania bed ze schematem primerow i tworzeniem na jego podstawie slownika. Offset determinuje
    o ile bp rozszerzamy amplikon w strone 5' i 3'. Z default offset to 0 i raczej tego nie bede zmienial
    :param bed: sciezka do pliku bed z primerami. Opisane w dokumentacji
    :param bed_offset: int o ile rozszerzamy amplikon. Default 1 zakladamym ze uzytkownik poda poprawnego bed-a a
    my sami rozszerzamy go o 1
    :return: 5  slownikow. Slowniki slownik_amplikonow_with_alt_outer i slownik_amplikonow_with_alt_inner maja jako
    klucz numer amplikonu (1,2,3...). Ktore same sa slownikami z tylko dwoma kluczami (LEFT i RIGHT). Wartosciami tych
    podslownikow jest lista. Najczesciej jedno elementowa zawierajaca granice primeru (zewnetrzna lub wewnetrzna w
    zaleznosci od slownika). Jest lista poniewaz trzymamy tam informacje o ewentualnych primerach alt. tak wiec
    {1}:{LEFT}:[1,3] oznacza ze w slowniku jest amplikon 1 i jego primer left zaczyna sie (lub konczy) na pozycji 1,
    ale jest tez primer alternatywny z poczatkiem (koncem na pozycji 3). Uwaga w slownilu zachowujemy sytuacje z bed-a
    w ktorym zakres jest polotwrty <)
    3 ostatnie slowniki to tylko rzeczy do statystyk gdzie jako klucz mamy numer amplikonu i wartosc zawsze 0
    na tym etapie

    Stara pomoc
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


    """
    slownik_amplikonow_with_alt_outer={}  # slowanik z zewnetrznymi granicami amplikonu (czyl 5' primer left i 3' right)
    slownik_amplikonow_with_alt_inner={}  # to samo ale granice wewnetrzne (3' left i 5' right)
    slownik_amplikonow_uzycie={}  # slownik z iloscia odczytow mapujacych sie wewnatrz danego amplikonu
    slownik_amplikonow_uzycie_left={}  # slownik z iloscia odczytow mapujacych sie na primer left danego amplikonu
    slownik_amplikonow_uzycie_right={}  # slownik z iloscia odczytow mapujacych sie na primer right danego amplikonu

    with open(bed) as f:
        for line in f:
            line = line.split()
            # pole 3 moze miec 3 albo 4 elementy bo moga wystepowac amplikony alt
            try:
                _, numer, kierunek, alt = line[3].split('_')
            except:
                _, numer, kierunek = line[3].split('_')

            numer = int(numer)

            if numer not in slownik_amplikonow_with_alt_outer.keys():

                slownik_amplikonow_with_alt_outer[numer] = {} # Ten slownik trzyma zewnetrzne granice amplikonow
                slownik_amplikonow_with_alt_inner[numer] = {} # Ten slownik trzyma wewnetrzne granice amplikonow

                slownik_amplikonow_with_alt_outer[numer]['LEFT'] = []
                slownik_amplikonow_with_alt_outer[numer]['RIGHT'] = []

                slownik_amplikonow_with_alt_inner[numer]['LEFT'] = []
                slownik_amplikonow_with_alt_inner[numer]['RIGHT'] = []

                slownik_amplikonow_uzycie[numer] = 0
                slownik_amplikonow_uzycie_left[numer] = 0
                slownik_amplikonow_uzycie_right[numer] = 0

            if 'LEFT' == kierunek :
                slownik_amplikonow_with_alt_outer[numer]['LEFT'].append(int(line[1]) - bed_offset)
                slownik_amplikonow_with_alt_inner[numer]['LEFT'].append(int(line[2]))
            elif 'RIGHT' == kierunek:
                slownik_amplikonow_with_alt_outer[numer]['RIGHT'].append(int(line[2]) + bed_offset)
                slownik_amplikonow_with_alt_inner[numer]['RIGHT'].append(int(line[1]))
            else:
                print('Nie rozpoznany kierunek')

    return slownik_amplikonow_with_alt_outer, slownik_amplikonow_with_alt_inner, \
        slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right, line[0]

def get_amplikon_coverage_in_windows(slownik_amplikonow, slownik_cap, szerokosc_okna=100):
    """
    Funkcja tworzy slownik z pokryciami w oknach. Wartosc pokrycia wynoi tyle co stala wartosc pobrana ze slownika
    pokrycia
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
                slownik_pokrycia_local[klucz][f"{initial_split[i]}_{initial_split[i + 1] - 1}"] = slownik_cap[klucz]
            except IndexError:
                # jestesmy na ostatnim elemencie listy
                # nie ma wiec wartosci o indekscie i+1
                # usuwamy poprzeni wpis i zamieniamy go tak by konczyl sie na wartosci end
                del(slownik_pokrycia_local[klucz][f"{initial_split[i-1]}_{initial_split[i] - 1}"])
                slownik_pokrycia_local[klucz][f"{initial_split[i-1]}_{end}"] = slownik_cap[klucz]

    return slownik_pokrycia_local

def filter_reads(initial_bam, final_bam, dlugosci_segmentow, statystyki, mapq=30, min_overlap=0.6,
                 min_alignment_overlap=0.65):
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
    :param min_overlap: jaka minimalna dlugosc musi miec odczyt (wyrazany jako procent dlugosci segmentu na
    jaki sie mapuje)

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
        # przekazemy po prostu dlugosc pierwszego amplikonu
        dlugosc_segmentu = dlugosci_segmentow


        if read.mapq < mapq:
            done = False
            statystyki.write(f"{read.qname}\tPoor mapq\n")
        elif (reference_end - reference_start) < (dlugosc_segmentu * min_alignment_overlap):
            done = False
            statystyki.write(f"{read.qname}\tAligned region is too short\n")
        elif read.rlen < dlugosc_segmentu * min_overlap:
            done = False
            statystyki.write(f"{read.qname}\tRead is too long {read.rlen}, segment has {dlugosc_segmentu}\n")


        if done:
            statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tPassed QC\n")
            all_reads_pass.write(read)

    all_reads_pass.close()
    pysam.collate('-o', final_bam, 'tmp.bam')
    os.remove('tmp.bam')

    return True

def write_reads_strict_inner(initial_bam, final_bam, reject_bam,  statystyki, slownik_amplikonow_with_alt_outer,
                      slownik_amplikonow_with_alt_inner, slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left,
                      slownik_amplikonow_uzycie_right, cap = 1000):
    """
    Ta funkcja wyciaga odczyty ktore zaczynaja sie i koncza w w primer. Dodatkowo skoro ta funkcja jest zawsze wykonywana
    na tym etapie wyrzucamy read o jakosci mapowania ponizej 30 i dlugosci ponizej length. To ma na celu przyspieszenie
    obliczen bo takie ready i tak nigdzie nie trafia
    ----------PRIMERL--------------------------PRIMERR----------
    ------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx-------------
    :param initial_bam:
    :param final_bam:
    :param reject_bam:
    :param statystyki:
    :param slownik_amplikonow_with_alt_outer:
    :param slownik_amplikonow_with_alt_inner:
    :param slownik_amplikonow_uzycie:
    :param slownik_amplikonow_uzycie_left:
    :param slownik_amplikonow_uzycie_right:
    :param cap:
    :return:
    """

    all_reads = pysam.AlignmentFile(initial_bam, "rb", require_index=False)
    pass_reads = pysam.AlignmentFile(final_bam, "wb", template=all_reads)
    reject_reads = pysam.AlignmentFile(reject_bam, "wb", template=all_reads)
    # uwaga pysam bierze koordynaty z bam-a z 1-index na 0-index


    for read in all_reads:
        done = False
        reference_start = read.reference_start
        reference_end = read.reference_end

        for klucz in slownik_amplikonow_with_alt_outer:
            if done:
                continue
            for indeks_start in range(len(slownik_amplikonow_with_alt_outer[klucz]['LEFT'])):
                for indeks_end in range(len(slownik_amplikonow_with_alt_outer[klucz]['RIGHT'])):
                    begin_amplikon_L = slownik_amplikonow_with_alt_outer[klucz]['LEFT'][indeks_start]
                    end_amplikon_L = slownik_amplikonow_with_alt_inner[klucz]['LEFT'][indeks_start]

                    begin_amplikon_R = slownik_amplikonow_with_alt_inner[klucz]['RIGHT'][indeks_end]
                    end_amplikon_R = slownik_amplikonow_with_alt_outer[klucz]['RIGHT'][indeks_end]

                    if begin_amplikon_L <= reference_start < end_amplikon_L \
                            and end_amplikon_R > reference_end >= begin_amplikon_R\
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] <= cap:
                        slownik_amplikonow_uzycie_left[klucz] += 1
                        slownik_amplikonow_uzycie_right[klucz] += 1
                        slownik_amplikonow_uzycie[klucz] += 1
                        pass_reads.write(read)
                        done = True
                        statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon "
                                         f"{klucz } both primers \n")
                    elif begin_amplikon_L <= reference_start < end_amplikon_L \
                            and end_amplikon_R > reference_end >= begin_amplikon_R\
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] > cap:
                        statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon "
                                         f"{klucz } both primers above cap \n")
                        done = True
                        # Dodaje tez info o rzeczywistych uzyciach nawet jesli bam uzywany do analizy ma mniejszy cap
                        slownik_amplikonow_uzycie_left[klucz] += 1
                        slownik_amplikonow_uzycie_right[klucz] += 1
                        slownik_amplikonow_uzycie[klucz] += 1

        if not done:
            reject_reads.write(read)

    pass_reads.close()
    reject_reads.close()

    name=final_bam.split('.')[0]
    pysam.sort('-o', f'{name}_sort.bam', final_bam)
    pysam.index(f'{name}_sort.bam')


    return slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right

def write_reads_overshot(initial_bam, final_bam, reject_bam, statystyki, slownik_amplikonow_with_alt_outer,
                         slownik_amplikonow_with_alt_inner, slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left,
                         slownik_amplikonow_uzycie_right, overshoot,  cap):
    """
    Ta funkcja zapisuje ready ktore przestrzelowne sa o 10 w stosunku do pozycji amplikonu (jednego badz obu)
    drugiego primera
    ----------PRIMERL--------------------------PRIMERR--------------
    --------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx----------- (przestrzelone oba primery na zewnatrz)
    --------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx----------------- (przestrzelony primer L na zewnatrz)
    -------------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx--------- (przestrzelony primer R na zewnatrz)
    -------------------xxxxxxxxxxxxxxxxxxxxx------------------------ (przestrzelone  L i R wewnatrz o wartosc overshot)
    -------------------------XXXXXXXXX------------------------------ (TE NIE WCHODZA BO SA ZA DALEKO OD PRIMERU)
    Granica 10 wynika jest empiryczna ale mozna sie podpudowac ze jest to polowa odlegosci do pozycji ktora jest miedzy
    primerami sasiednich amplikonow (primery roznych amplikonow sa oddalone o ok 40 bp)
    :param initial_bam:
    :param final_bam:
    :param statystyki:
    :param slownik_amplikonow_with_alt_outer:
    :param slownik_amplikonow_with_alt_inner:
    :param slownik_amplikonow_uzycie:
    :param slownik_amplikonow_uzycie_left:
    :param slownik_amplikonow_uzycie_right:
    :param cap:
    :return:
    """
    all_reads = pysam.AlignmentFile(initial_bam, "rb",  require_index=False)
    pass_reads = pysam.AlignmentFile(final_bam, "wb", template=all_reads)
    reject_reads = pysam.AlignmentFile(reject_bam, "wb", template=all_reads)
    # uwaga pysam bierze koordynaty z bam-a z 1-index na 0-index

    for read in all_reads:
        done = False
        reference_start = read.reference_start
        reference_end = read.reference_end

        for klucz in slownik_amplikonow_with_alt_outer:
            if done:
                continue
            for indeks_start in range(len(slownik_amplikonow_with_alt_outer[klucz]['LEFT'])):
                for indeks_end in range(len(slownik_amplikonow_with_alt_outer[klucz]['RIGHT'])):
                    begin_amplikon_L = slownik_amplikonow_with_alt_outer[klucz]['LEFT'][indeks_start]
                    end_amplikon_L = slownik_amplikonow_with_alt_inner[klucz]['LEFT'][indeks_start]

                    begin_amplikon_R = slownik_amplikonow_with_alt_inner[klucz]['RIGHT'][indeks_end]
                    end_amplikon_R = slownik_amplikonow_with_alt_outer[klucz]['RIGHT'][indeks_end]

                    if ( begin_amplikon_L - overshoot ) <= reference_start < (end_amplikon_L + overshoot)  \
                            and (end_amplikon_R + overshoot) > reference_end >= (begin_amplikon_R - overshoot) \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] <= cap:

                        slownik_amplikonow_uzycie_left[klucz] += 1
                        slownik_amplikonow_uzycie_right[klucz] += 1
                        slownik_amplikonow_uzycie[klucz] += 1

                        pass_reads.write(read)
                        done = True
                        statystyki.write(
                            f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz} overshoot\n")
                    elif ( begin_amplikon_L - overshoot ) <= reference_start < (end_amplikon_L + overshoot)  \
                            and (end_amplikon_R + overshoot)  > reference_end >= (begin_amplikon_R - overshoot) \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] > cap:
                        done = True
                        statystyki.write(
                            f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz} overshoot  cap\n")

        if not done:
            reject_reads.write(read)

    pass_reads.close()
    reject_reads.close()

    name = final_bam.split('.')[0]
    pysam.sort('-o', f'{name}_sort.bam', final_bam)
    pysam.index(f'{name}_sort.bam')

    return slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right

def write_reads_fusion_strict(initial_bam, final_bam, statystyki, primer_left_outer, primer_left_inner,
                              primer_middle_left_outer,  primer_middle_left_inner, primer_middle_right_outer,
                              primer_middle_right_inner, primer_right_outer, primer_right_inner, used, case,
                              slownik_amplikonow_uzycie, klucz,  cap, overshoot):
    """
    To jest modyfikacja funkcji write_reads_fusion w ktorej zakladamy ze smieci musza zaczynac sie w primerach
    tzn albo jest fuzja -2 i 0 albo 0 + 2 albo wrecz -2 do + 2 i kazda z tych mozliwosci rozwazamy, dla amplikonow pierwsz
    ego i ostatniego oczywiscie czesc z tych mozliwosci jest zablokowana. Case 1 to przypadek fdzie nie ma primera -2,
    Case 2 gdzie nie ma primera +2 , case 3 gdzie jest primer -2 i +2 . Zeby bylo juz latwiej gdy mamy alty w primerach
    to do tej funkcji dajemy po prosy maksymalne zasiegi danego primera czyli zwykly primer + alt i zmienne primer_left_outer itd
    sa po prostu int-ami
    :param initial_bam:
    :param final_bam:
    :param statystyki:
    :param primer_left_outer:
    :param primer_left_inner:
    :param primer_middle_left_outer:
    :param primer_middle_left_inner:
    :param primer_middle_right_outer:
    :param primer_middle_right_inner:
    :param primer_right_outer:
    :param primer_right_inner:
    :param used:
    :param case:
    :param slownik_amplikonow_uzycie:
    :param klucz:
    :param cap:
    :param overshoot:
    :return:
    """
    my_bam = pysam.AlignmentFile(initial_bam, "rb", require_index=False)
    my_bam_out = pysam.AlignmentFile(final_bam, "wb", template=my_bam)

    amplikon_ilosc = slownik_amplikonow_uzycie[klucz]

    # dla skrajnych amplikonow tych kluczy nie bedzie wiec zrobimy je dummy
    try:
        amplikon_ilosc_left_left = slownik_amplikonow_uzycie[klucz - 2]
    except KeyError:
        amplikon_ilosc_left_left = 0
    try:
        amplikon_ilosc_left = slownik_amplikonow_uzycie[klucz - 1]
    except KeyError:
        amplikon_ilosc_left = 0

    try:
        amplikon_ilosc_right = slownik_amplikonow_uzycie[klucz + 1]
    except KeyError:
        amplikon_ilosc_right = 0

    try:
        amplikon_ilosc_right_right = slownik_amplikonow_uzycie[klucz + 2]
    except KeyError:
        print(f'W slowniku nie ma klucza {klucz + 2}')
        amplikon_ilosc_right_right = 0

    # skoto obejmuje cale amplikony to podbijam tez ich uzycie
    with open(statystyki, 'w') as f:
        # lecimy po parach odczytów
        for read in my_bam:
            if read.qname in used.keys():
                #odczyt zostal juz wykorzystany pomijam go wiec
                continue

            reference_start = read.reference_start
            reference_end = read.reference_end
            # w odroznieniu od funkcji wyzej tu nie dajemy slownika amplikonow
            # a wybrane listy z zakresami


            if (case == 1 or case == 3) and (primer_middle_left_outer - overshoot ) <= reference_start < (primer_middle_left_inner + overshoot) \
                    and (primer_right_inner - overshoot ) <= reference_end < ( primer_right_outer + overshoot ) and amplikon_ilosc <= cap:
                    # read jest od lewego primera swojego do prawego priemera +2
                my_bam_out.write(read)
                used[read.qname] = ''
                amplikon_ilosc += 1
                amplikon_ilosc_left += 1
                amplikon_ilosc_left_left +=1
                f.write(
                    f"{read.qname}\t{read.reference_start}\t{read.reference_end}\ttwo_amplicons_{klucz}\n")
                continue
            elif (case == 2 or case == 3) and (primer_left_outer - overshoot ) <= reference_start < (primer_left_inner + overshoot ) and \
                    (primer_middle_right_inner - overshoot ) <= reference_end < (primer_middle_right_outer + overshoot ) and amplikon_ilosc <= cap:
                # read jest od lewego primera -2 do swojego primer prawgo
                my_bam_out.write(read)
                used[read.qname] = ''
                amplikon_ilosc += 1
                amplikon_ilosc_right += 1
                amplikon_ilosc_right_right += 1
                f.write(
                    f"{read.qname}\t{read.reference_start}\t{read.reference_end}\ttwo_amplicons_{klucz}\n")
                continue
            elif case == 3 and (primer_left_outer - overshoot) <= reference_start < (primer_left_inner + overshoot ) and \
                    (primer_right_inner - overshoot )<= reference_end < (primer_right_outer + overshoot ) and amplikon_ilosc <= cap:
                # odczyt obejmuje region od amplikonu -2 do +2 a wiec bardzo dlugi odcinek
                # bo "ominieto" oba primery wypadnietego amplikonu
                my_bam_out.write(read)
                used[read.qname] = ''
                amplikon_ilosc += 1
                amplikon_ilosc_left += 1
                amplikon_ilosc_left_left += 1
                amplikon_ilosc_right += 1
                amplikon_ilosc_right_right += 1
                f.write(
                    f"{read.qname}\t{read.reference_start}\t{read.reference_end}\ttwo_amplicons_{klucz}\n")
                continue


    my_bam_out.close()

    #final_bam_sort = final_bam.split('.')[0]
    #pysam.sort("-o", f"{final_bam_sort}_sort.bam", final_bam)
    #pysam.index(f"{final_bam_sort}_sort.bam")

    # updatujemy slownik uzycia o obecne w nim oryginalne klucze
    if ( klucz - 2) in slownik_amplikonow_uzycie.keys():
        slownik_amplikonow_uzycie[klucz - 2 ] = amplikon_ilosc_left_left
    if (klucz - 1) in slownik_amplikonow_uzycie.keys():
        slownik_amplikonow_uzycie[klucz - 1] = amplikon_ilosc_left
    if ( klucz) in slownik_amplikonow_uzycie.keys():
        slownik_amplikonow_uzycie[klucz] = amplikon_ilosc
    if (klucz + 1) in slownik_amplikonow_uzycie.keys():
        slownik_amplikonow_uzycie[klucz + 1] = amplikon_ilosc_right
    if (klucz + 2) in slownik_amplikonow_uzycie.keys():
        slownik_amplikonow_uzycie[klucz + 2] = amplikon_ilosc_right_right

    return slownik_amplikonow_uzycie, used

def write_reads_partstrict_inner(initial_bam, final_bam, reject_bam, statystyki, slownik_amplikonow_with_alt_outer,
                      slownik_amplikonow_with_alt_inner, slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left,
                      slownik_amplikonow_uzycie_right, cap = 1000, overshoot = 0, usage = 0.6):
    """
    Ta funkcja wyciaga odczyty ktorych mapowanie zaczyna w jednym z primerow. Mniej cenne niz write_reads_strict_inner
    ale tez dosc dobrze wiadomo skad pochodza odczyty
    Przywrocilem ta funkcje bo w sampu 9 wyraznie korzystaja z takich readow, gorzej ze nie sa zbyt konsystentni
    DAjemy empitycznie ze read obejmuje jednak co najmniej polowe amplikonu
    ----------PRIMERL--------------------------PRIMERR----------
    ------------xxxxxxxxxxxxxxxxxxxxxxxxxx----------------------
    --------------------xxxxxxxxxxxxxxxxxxxxxxxxxx--------------
    ALA takie odczyty NIE MOGA
    ----------PRIMERL------------------PRIMERLinnyklucz--------------PRIMERR
    ------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx------------------------------
    ALE moga byc tak
    --------------xxxxxxxxxxxxxxxxx-----------------------------------------
    po warunkiem ze to ponad 0.5 uzycia maplikonu. Ta funkcja jest strasznie pod EQA

    zrobic czego takiego
    Ba tym etapie filtruje tez ready ktore sa po prosry za krotkiw
    :param initial_bam:
    :param final_bam:
    :param reject_bam:
    :param statystyki:
    :param slownik_amplikonow_with_alt_outer:
    :param slownik_amplikonow_with_alt_inner:
    :param slownik_amplikonow_uzycie:
    :param slownik_amplikonow_uzycie_left:
    :param slownik_amplikonow_uzycie_right:
    :param cap:
    :return:
    """
    all_reads = pysam.AlignmentFile(initial_bam, "rb", require_index=False)
    pass_reads = pysam.AlignmentFile(final_bam, "wb", template=all_reads)
    reject_reads = pysam.AlignmentFile(reject_bam, "wb", template=all_reads)
    # uwaga pysam bierze koordynaty z bam-a z 1-index na 0-index


    for read in all_reads:
        done = False
        reference_start = read.reference_start
        reference_end = read.reference_end

        for klucz in slownik_amplikonow_uzycie:
            if done:
                continue
            for indeks_start in range(len(slownik_amplikonow_with_alt_outer[klucz]['LEFT'])):
                for indeks_end in range(len(slownik_amplikonow_with_alt_outer[klucz]['RIGHT'])):
                    begin_amplikon_L = slownik_amplikonow_with_alt_outer[klucz]['LEFT'][indeks_start]
                    end_amplikon_L = slownik_amplikonow_with_alt_inner[klucz]['LEFT'][indeks_start]

                    begin_amplikon_R = slownik_amplikonow_with_alt_inner[klucz]['RIGHT'][indeks_end]
                    end_amplikon_R = slownik_amplikonow_with_alt_outer[klucz]['RIGHT'][indeks_end]

                    # jesli read zaczyna sie w lewym primer to 'musi' przejsc lewy primer kolejnego amplikonu + offset

                    try:
                        end_amplikon_L_plus = slownik_amplikonow_with_alt_inner[klucz + 1]['LEFT'][0]
                        begin_amplikon_L_plus = slownik_amplikonow_with_alt_outer[klucz + 1]['LEFT'][0]
                    except:
                        end_amplikon_L_plus = 0
                        begin_amplikon_L_plus = 0

                    # jesli read zaczyna sie w primer prawy to musi przejsc prawy primer amplikonu -1 - offset
                    try:
                        begin_amplikon_R_minus =  slownik_amplikonow_with_alt_inner[klucz - 1 ]['RIGHT'][0]
                        end_amplikon_R_minus = slownik_amplikonow_with_alt_outer[klucz - 1]['RIGHT'][0]
                    except:
                        begin_amplikon_R_minus = 10000000
                        end_amplikon_R_minus = 10000000

                    # Orygianly warunek w EQA2024 byl taki
                    # ze odczyt zaczynajcy sie w primer musi siegac poza primer nastepnego amplikonu
                    # z zachowaniem kierunku to znaczy
                    # jesli startuje z primer'u lewego musze siegnac poza kolejny lewy primer (czyli amplikon +1)
                    # jesli startuje z primeru prawego musze siegnac poza kolejny prawy primer (jako ze sie
                    # "cofam" to bedzie to primer amplikonu -1

                    # Nie wiem skad ten warunek byl ale jest to de facto constrain na dlugosc
                    # Mozliwe tez ze samtoools ampliconclip sie wtedy gubi i cos zle przycina

                    #usune ten warunk jednak ale z zachowaniem ze dla artica taki odczyt za krotki jest niepozadany
                    # stad uzycie bedzie bardzo niskie

                    if (begin_amplikon_L - overshoot) <= reference_start < (end_amplikon_L + overshoot )\
                            and begin_amplikon_R >= reference_end \
                            and ( reference_end > end_amplikon_L_plus or  reference_end < begin_amplikon_L_plus) \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] <= cap:

                        amplicon_zakres = set(range(begin_amplikon_L, end_amplikon_R))
                        read_zakres = set(range(reference_start,reference_end))
                        common_zakres = len(read_zakres.intersection(amplicon_zakres)) / float(len(amplicon_zakres))
                        if common_zakres >= usage:
                            slownik_amplikonow_uzycie_left[klucz] += 1
                            slownik_amplikonow_uzycie[klucz] += 1
                            pass_reads.write(read)
                            done = True
                            statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz} only left primer\n")
                        else:
                            done = True
                    elif end_amplikon_L < reference_start \
                            and (end_amplikon_R  + overshoot ) > reference_end >= (begin_amplikon_R - overshoot )\
                            and ( reference_start < begin_amplikon_R_minus or  reference_start > end_amplikon_R_minus) \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] <= cap:
                        amplicon_zakres = set(range(begin_amplikon_L, end_amplikon_R))
                        read_zakres = set(range(reference_start, reference_end))
                        common_zakres = len(read_zakres.intersection(amplicon_zakres)) / float(len(amplicon_zakres))
                        if common_zakres >= usage:
                            slownik_amplikonow_uzycie_right[klucz] += 1
                            slownik_amplikonow_uzycie[klucz] += 1
                            pass_reads.write(read)
                            done = True
                            statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz} only right primer\n")
                        else:
                            done = True
        if not done:
            reject_reads.write(read)

    pass_reads.close()
    reject_reads.close()

    final_bam_sort = final_bam.split('.')[0]
    pysam.sort("-o", f"{final_bam_sort}_sorted.bam", final_bam)
    pysam.index(f"{final_bam_sort}_sorted.bam")

    return slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right


def write_reads_midnight(initial_bam, final_bam, reject_bam, statystyki, slownik_amplikonow_with_alt_outer,
                      slownik_amplikonow_uzycie, cap = 1000):
    """
    Funkcja ktora bierze smiecoi z midnight, tzn odczyty ktores sa wewnatrz amplikonow, ale nie obejmuja primerow,
    to sa odczyty po filtrowaniu na dlugosc
    :param initial_bam:
    :param final_bam:
    :param reject_bam:
    :param statystyki:
    :param slownik_amplikonow_with_alt_outer:
    :param slownik_amplikonow_uzycie:
    :param cap:
    :return:
    """

    all_reads = pysam.AlignmentFile(initial_bam, "rb", require_index=False)
    pass_reads = pysam.AlignmentFile(final_bam, "wb", template=all_reads)
    reject_reads = pysam.AlignmentFile(reject_bam, "wb", template=all_reads)
    # uwaga pysam bierze koordynaty z bam-a z 1-index na 0-index


    for read in all_reads:
        done = False
        reference_start = read.reference_start
        reference_end = read.reference_end

        for klucz in slownik_amplikonow_with_alt_outer:
            if done:
                continue
            for indeks_start in range(len(slownik_amplikonow_with_alt_outer[klucz]['LEFT'])):
                for indeks_end in range(len(slownik_amplikonow_with_alt_outer[klucz]['RIGHT'])):
                    begin_amplikon_L = slownik_amplikonow_with_alt_outer[klucz]['LEFT'][indeks_start]
                    end_amplikon_R = slownik_amplikonow_with_alt_outer[klucz]['RIGHT'][indeks_end]

                    if begin_amplikon_L  <= reference_start  \
                            and end_amplikon_R  > reference_end \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] <= cap:
                        slownik_amplikonow_uzycie[klucz] += 1
                        pass_reads.write(read)
                        done = True
                        statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz } midnight case \n")
                    elif begin_amplikon_L  <= reference_start  \
                            and end_amplikon_R  > reference_end \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] > cap:
                        done = True
                        statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz } midnight case cap \n")


        if not done:
        # probowalem read nie trafil do wynikowego skryptu
        # dajemy mu 2-ga szanse w kolejnym pass
            reject_reads.write(read)

    pass_reads.close()
    reject_reads.close()

    final_bam_sort = final_bam.split('.')[0]
    pysam.sort("-o", f"{final_bam_sort}_sort.bam", final_bam)
    pysam.index(f"{final_bam_sort}_sort.bam")
    return slownik_amplikonow_uzycie


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

def update_slownik_pokrycia(slownik_pokrycia, reference_start, reference_end):
    """
    Funkcja do updatu wartosci w slowniku z pokryciem.
    :param slownik_pokrycia: slownik w ktorym kluczami sa zakresy okien (np '100_200') a jedyna wartoscia int z
     aktualnym pokryciem. To wycinek slownika tworzonego przez funkcje get_empty_amplikon_coverage_in_windows.
    :param reference_start: Poczatek mapowania pierwszego z odczytow z pary
    :param reference_end: Koniec mapowania pierwszego z odczytow z pary
    :return: update'owany slownik z pokryciem o strukturze identycznej jak ten podawany jako input
    """

    odczyt_zakres = set(range(reference_start, reference_end))

    for klucz in slownik_pokrycia.keys():
        # iterujemy po slowniku pokryciu klucze to zakresy np. '0_99', '100-199 itd ...

        klucz_start = int(klucz.split('_')[0])
        klucz_end = int(klucz.split('_')[1])
        okno_zakres = set(range(klucz_start, klucz_end))

        #jaka frakcja okna jest obejmowana przez odczyt i pare tego odczytu
        pokrycie_read = len(okno_zakres.intersection(odczyt_zakres)) / (klucz_end - klucz_start)

        #update wartosci slownika
        slownik_pokrycia[klucz] += pokrycie_read

    return slownik_pokrycia


def get_inneramplicon_reads_window(reads, reads_pass, reads_reject, slownik_pokrycia,
                            slownik_amplikonow_outer, slownik_amplikonow_inner, uzycie_left, uzycie_right,
                            cap, statystyki):
    """
    Funkcja do wyrownania pokrycia tylko odczytami mapujacymi sie na jeden amplikon
    :param reads: Sciezka do pliku ze wszystkimi odczytami
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

    all_reads_forward = pysam.AlignmentFile(reads, 'rb', require_index=False)


    pass_reads = pysam.AlignmentFile(f'{reads_pass}', 'wb', template=all_reads_forward)
    reject_reads = pysam.AlignmentFile(f'{reads_reject}','wb', template=all_reads_forward)


    for odczyt in all_reads_forward:
        # gdzie mapuje sie odczyt i jego mate
        reference_start = odczyt.reference_start
        reference_end = odczyt.reference_end
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

            if reference_start >= amplikon_start and reference_end < amplikon_end:
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
                    if amplikon_start_inner > reference_start >= amplikon_start:
                        # odczyt zaczyna sie w lewym primer
                        uzycie_left[amplikon] += 1
                    if amplikon_end_inner <= reference_end < amplikon_end:
                        # odczyt konczy sie w amplikonie prawym
                        uzycie_right[amplikon] += 1

                    # update pokrycia tym odczytem
                    slownik_pokrycia[amplikon] = update_slownik_pokrycia(slownik_pokrycia[amplikon], reference_start,
                                                                         reference_end)

                    # zapisanie info do pliku
                    statystyki.write(f'Odczyt {odczyt.qname} z zakresow {reference_start} {reference_end} '
                                     f'inside amplikonu {amplikon} below cap\n')
                else:
                    # odczyt jest nieprzydatny zapisujemy info do pliku
                    statystyki.write(f'Odczyt {odczyt.qname} z zakresow {reference_start} {reference_end} '
                                     f' inside amplikonu {amplikon} above cap\n')

        # zbadalem wszystkie amplikony dla danego odczytu zapisuje wynik do plikow
        if done:
            # zapisanie pary odczytow do pliku first pass
            pass_reads.write(odczyt)
        if not one_amplicon and not done:
            # odczyt jest poza pojednyczym amplikonem i nie trafil co oczywiste do pliku first_pass

            statystyki.write(f'Odczyt {odczyt.qname} z zakresow {reference_start} {reference_end} '
                             f'not in amplicon\n')
            reject_reads.write(odczyt)


    # zamkniecie plikow i ich sortowanie z indeksowaniem
    pass_reads.close()
    reject_reads.close()

    pass_sorted_name = f'{reads_pass.split(".")[0]}_sort.bam'
    reject_sorted_name = f'{reads_reject.split(".")[0]}_sort.bam'

    pysam.sort("-o", pass_sorted_name, reads_pass)
    pysam.index(pass_sorted_name)

    pysam.sort("-o", reject_sorted_name, reads_reject)
    pysam.index(reject_sorted_name)
    return slownik_pokrycia, uzycie_left, uzycie_right, pass_sorted_name, reject_sorted_name

if __name__ == '__main__':
    all_read = sys.argv[1]  # posortowany output z minimap-a
    bed = sys.argv[2]  # plik bed z primerami
    bed_offset = int(sys.argv[3])  # rozszerzanie amplikonu w kierunku 5' i 3' na etapie strict
    cap = int(sys.argv[4])  # cap na ILOSC odczytow mapujacych sie na konkretny amplikon
    length_min = float(sys.argv[5])  # minmalna dlugosc odczytu wyrazona jako frakcja dlugosci amplikonu
    mapq = int(sys.argv[6])
    extra_bed_offset = int(sys.argv[7])
    cap_smieci=int(sys.argv[8])  # cap na smieci czyli ready nie obejmujace dwoch primerow
    statystyki = open('Statystyki.txt', 'w')

    #1 Slownik z ampikonami i uzyciami amplkionow
    slownik_amplikonow_with_alt_outer, slownik_amplikonow_with_alt_inner, \
        slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, \
        slownik_amplikonow_uzycie_right, ref_name = read_amplicon_scheme(bed=bed, bed_offset=bed_offset)



    midnight = False # midnigh ma tagmentacje i odczyty nie beda zawsze obejmowaly calego amplikonu
    dlugosc_pierwszego_amplikonu = (slownik_amplikonow_with_alt_outer[1]['RIGHT'][0] -
                                    slownik_amplikonow_with_alt_outer[1]['LEFT'][0])

    # podmiana cap smiecie na cap w przypadku dlugich amplikonow z midnight, ktore i
    # tak sa ciete na krotsze odcinki przez tagmentaze
    if (dlugosc_pierwszego_amplikonu) > 600:
        cap_smieci = cap
        midnight = True
        length_min = 0.2  # pozwalamy aby odczyty mialy tylko 20% dlugosci amplikonu
        # bo w protokole Tomka jest duzo fragmentow


    # 1 standardowe filtrowanie
    filter_reads(initial_bam=all_read,
                 final_bam='all_reads_clean.bam',
                 dlugosci_segmentow=dlugosc_pierwszego_amplikonu,
                 statystyki=statystyki,
                 min_overlap=length_min,
                 min_alignment_overlap=0.49,
                 mapq=mapq)

    # Na tym etapie odczyty sa juz zrandomizowane !
    #2 Pierwszy pass odczyty strict
    # Na tym etapie okna sa nieistotne bo odczyt ma objac caly amplikon wiec cap bedzie rowny
    # sredniemu pokryciu
    slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, \
        slownik_amplikonow_uzycie_right = write_reads_strict_inner(initial_bam='all_reads_clean.bam',
                                                                   final_bam='reads_inner_strict.bam',
                                                                   reject_bam='reject_first_pass.bam',
                                                                   statystyki=statystyki,
                                                                   slownik_amplikonow_with_alt_outer=slownik_amplikonow_with_alt_outer,
                                                                   slownik_amplikonow_with_alt_inner=slownik_amplikonow_with_alt_inner,
                                                                   slownik_amplikonow_uzycie=slownik_amplikonow_uzycie,
                                                                   slownik_amplikonow_uzycie_left=slownik_amplikonow_uzycie_left,
                                                                   slownik_amplikonow_uzycie_right=slownik_amplikonow_uzycie_right,
                                                                   cap=cap)

    # pysam.index('reads_inner_strict.bam')

    # pysam.sort('-o', 'reject_first_pass_sort.bam', 'reject_first_pass.bam')
    # pysam.index("reject_first_pass_sort.bam")

    ### PO etapie ###
    ### Mamy posortowany i zindekwoany plik z odczytami obejmujacymi pelen amplikon
    ### Oraz dalej nieposortownay i niezindeksowany reject_first_pass.bam
    ###  ####

    #3 Ready ktore przestrzeliy jeden lub oba primery beda oddzielnie trimowane z wieksza tolerancja przez ampliconclip
    # Te odczyty maja amplikon rozszerzony jeszcze bardziej o extra_bed_offset

    slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, \
        slownik_amplikonow_uzycie_right = write_reads_overshot(initial_bam='reject_first_pass.bam',
                                                               final_bam='reads_overshot.bam',
                                                               reject_bam='reject_second_pass.bam',
                                                               statystyki=statystyki,
                                                               slownik_amplikonow_with_alt_outer=slownik_amplikonow_with_alt_outer,
                                                               slownik_amplikonow_with_alt_inner=slownik_amplikonow_with_alt_inner,
                                                               slownik_amplikonow_uzycie=slownik_amplikonow_uzycie,
                                                               slownik_amplikonow_uzycie_left=slownik_amplikonow_uzycie_left,
                                                               slownik_amplikonow_uzycie_right=slownik_amplikonow_uzycie_right,
                                                               overshoot=extra_bed_offset,
                                                               cap=cap)


    ### PO etapie ###
    ### Mamy posortowany i zindekwoany plik z odczytami obejmujacymi pelen amplikon z tolerancja +/- ekstra bed offset
    ### Oraz dalej nieposortownay i niezindeksowany reject_second_pass.bam
    ###  ####



    #4 Ready smieci nie wiadomo z jakiego sa maplikonu ale mapuja sie w duzej mierze na amplikon o niskim uzyciu
    # moze pochodzi z fuzji sasiednich ampikonow gdy wypada primer (z powodu delecji, albo slabej hybrydyzacji)
    # albo bezposrednio z materialu wyjsciowego generalnie smietnik uzywany do podbicia coverage

    used = {}  # aby zapobiec dublowaniu sie odczytu gdy analizujemy N-razy input tworzymy slownik gdzie kluczami
    # sa uzyte odczyty
    # zeby pomoc amplikonom 1 i 2 oraz ostatniemu i przedostatniemu (ktore nie maja fuzji z ampikonami z obu stron
    #  modyfukjemy slownik_amplikonow_with_alt_outer
    # i slownik_amplikonow_with_alt_inner aby zawrieral primery dla -1 i -2 i + 1 i +2 takze dla nich
    # ale te symulowane granice sa identyczne ze skrajnymi primerami

    min_amplikon = min(slownik_amplikonow_with_alt_outer.keys()) # minimalna wartosc klczuami sa numery amplikonow  1,2,3
    max_amplikon = max(slownik_amplikonow_with_alt_outer.keys())

    slownik_amplikonow_with_alt_outer[min_amplikon - 1] =  slownik_amplikonow_with_alt_outer[min_amplikon]
    slownik_amplikonow_with_alt_inner[min_amplikon - 1] =  slownik_amplikonow_with_alt_inner[min_amplikon]

    slownik_amplikonow_with_alt_outer[min_amplikon - 2] =  slownik_amplikonow_with_alt_outer[min_amplikon]
    slownik_amplikonow_with_alt_inner[min_amplikon - 2] =  slownik_amplikonow_with_alt_inner[min_amplikon]

    slownik_amplikonow_with_alt_outer[max_amplikon + 1] =  slownik_amplikonow_with_alt_outer[max_amplikon]
    slownik_amplikonow_with_alt_inner[max_amplikon + 1] =  slownik_amplikonow_with_alt_inner[max_amplikon]

    slownik_amplikonow_with_alt_outer[max_amplikon + 2] =  slownik_amplikonow_with_alt_outer[max_amplikon]
    slownik_amplikonow_with_alt_inner[max_amplikon + 2] =  slownik_amplikonow_with_alt_inner[max_amplikon]

    # po tym may klucze -2, -1 oraz ostatni amplikon + 1 i ostatni amplikon +2

   #4 Odczyty powstale z fuzji sasiednich amplikonow z w granicach +/- EXTRA BED OFFSET
    for klucz, wartosc in slownik_amplikonow_uzycie.items():
        if wartosc < cap:
            # Dalej poruszamy sie w ramach odczytow obejmujacych pelen amplikon
            # Stad nie ma potrzeby wygladzania pokrycia w oknach

            # print(f'Analizuje fuzje ampliku {klucz} z sasiednimi amplikonami')

            initial_bam = "reject_second_pass.bam"
            final_bam = f'reads_two_amplicons_{klucz}.bam'
            statystyki_new = f'Statystyki_two_amplicons_{klucz}.txt'



            slownik_amplikonow_uzycie, used = write_reads_fusion_strict(initial_bam=initial_bam,
                                                                        final_bam=final_bam,
                                                                        statystyki=statystyki_new,
                                                                        primer_left_outer=min(slownik_amplikonow_with_alt_outer[(klucz - 2)]['LEFT']),
                                                                        primer_left_inner=max(slownik_amplikonow_with_alt_inner[(klucz - 2)]['LEFT']),
                                                                        primer_middle_left_outer=min(slownik_amplikonow_with_alt_outer[(klucz)]['LEFT']),
                                                                        primer_middle_left_inner=max(slownik_amplikonow_with_alt_inner[(klucz)]['LEFT']),
                                                                        primer_middle_right_outer=max(slownik_amplikonow_with_alt_outer[(klucz)]['RIGHT']),
                                                                        primer_middle_right_inner=min(slownik_amplikonow_with_alt_inner[(klucz)]['RIGHT']),
                                                                        primer_right_outer=max(slownik_amplikonow_with_alt_outer[(klucz + 2)]['RIGHT']),
                                                                        primer_right_inner=min(slownik_amplikonow_with_alt_inner[(klucz + 2)]['RIGHT']),
                                                                        slownik_amplikonow_uzycie=slownik_amplikonow_uzycie,
                                                                        klucz=klucz,
                                                                        used=used,
                                                                        case=3,
                                                                        overshoot = extra_bed_offset,
                                                                        cap=cap)

    ### PO etapie ###
    ### Mamy potenjalnie duzo plkow *two_amplicons* ktore NIE sa pososrtowane
    ### ale fragment shell-a jest zlaczy, pososrtuje i zamaskuje primery
    ###  ####

    #5 Odczyty ktore startuja w primer danego amplikonu, ale sa ZA krotkie
    # i nie dotarly do odpowiedniego primeru
    # Na tym etapie zachowujemy minimum przyzwoitosci i prefereujemy ready dluzsze nad krotszymi
    # krotkie beda na etapie ostatecznej desperacji

    # Odczyty te nie powinny podbijac pokrycia ponad cap_smieci

    # ale w tym wypadku chyba czas zaczac stosowac okna

    if midnight:
        usage = 0.3 #usage jest zbedny przy uzywaniu funkcji opartej o okna
    else:
        usage = 0.5
        cap = cap_smieci

    initial_bam = 'reject_second_pass.bam'
    odczyty_fragmenty = 'reads_partial_strict.bam'
    odczyty_fragmenty_odrzucone = 'reads_partial_strict_reject.bam'
    # Na tym etapie pokrycie w oknie powinno byc stale (bo odczyty obejmuja caly amplikon albo jego fuzne)
    # Tworzymy sobie wiec slownik pokrycia w oknach dopiero na tym etapie
    # i ewentualnie uzupelniamy
    # dla artic do cap_smieci (aby przekroczyc threshold na maskowanie)
    # dla midnight do "pelnej" wartosci
    del(slownik_amplikonow_with_alt_outer[min_amplikon - 1])
    del(slownik_amplikonow_with_alt_outer[min_amplikon - 2])
    del(slownik_amplikonow_with_alt_outer[max_amplikon + 1])
    del(slownik_amplikonow_with_alt_outer[max_amplikon + 2])
    slownik_pokrycia_w_oknach = get_amplikon_coverage_in_windows(slownik_amplikonow = slownik_amplikonow_with_alt_outer,
                                                                 slownik_cap = slownik_amplikonow_uzycie,
                                                                 szerokosc_okna = 50)
    print(slownik_pokrycia_w_oknach)
    slownik_pokrycia_w_oknach, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right, \
        first_pass_sorted_name, first_reject_sorted_name = \
        get_inneramplicon_reads_window(reads = initial_bam,
                                reads_pass=odczyty_fragmenty,
                                reads_reject=odczyty_fragmenty_odrzucone,
                                slownik_pokrycia=slownik_pokrycia_w_oknach,
                                slownik_amplikonow_outer=slownik_amplikonow_with_alt_outer,
                                slownik_amplikonow_inner=slownik_amplikonow_with_alt_inner,
                                uzycie_left=slownik_amplikonow_uzycie_left,
                                uzycie_right=slownik_amplikonow_uzycie_right,
                                cap=cap,
                                statystyki=statystyki)
    # print(slownik_pokrycia_w_oknach)
    # slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, \
    #     slownik_amplikonow_uzycie_right = write_reads_partstrict_inner(initial_bam=initial_bam,
    #                                                                    final_bam='reads_partial_strict.bam',
    #                                                                    reject_bam='reject_partial_strict_inner.bam',
    #                                                                    statystyki=statystyki,
    #                                                                    slownik_amplikonow_with_alt_outer=slownik_amplikonow_with_alt_outer,
    #                                                                    slownik_amplikonow_with_alt_inner=slownik_amplikonow_with_alt_inner,
    #                                                                    slownik_amplikonow_uzycie=slownik_amplikonow_uzycie,
    #                                                                    slownik_amplikonow_uzycie_left=slownik_amplikonow_uzycie_left,
    #                                                                    slownik_amplikonow_uzycie_right=slownik_amplikonow_uzycie_right,
    #                                                                    overshoot=extra_bed_offset,
    #                                                                    usage = usage,
    #                                                                    cap=cap_smieci)


    # To powinno byc zbedne, ale przetestowac na obserco
    # if midnight:
    #     slownik_amplikonow_uzycie =  write_reads_midnight(initial_bam = 'reject_partial_strict_inner_unsort.bam',
    #                                                       final_bam = 'reads_smieci.bam',
    #                                                       reject_bam = 'final_reject.bam',
    #                                                       statystyki = statystyki,
    #                                                       slownik_amplikonow_with_alt_outer = slownik_amplikonow_with_alt_outer,
    #                                                       slownik_amplikonow_uzycie = slownik_amplikonow_uzycie,
    #                                                       cap=cap_smieci)

    with open('Primer_usage.txt', 'w') as f:
        f.write('Segment\tPrimer_number\tPrimare_usage\n')
        for klucz, wartosc in slownik_amplikonow_uzycie.items():
            amplion_usage = slownik_amplikonow_uzycie_left[klucz] + slownik_amplikonow_uzycie_right[klucz]
            f.write(f'{ref_name}\t{klucz}\t{amplion_usage}\n')
:
    #pysam.index('reads_partial_strict.bam')






