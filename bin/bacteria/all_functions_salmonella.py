import os
import re
from Bio import SeqIO
import gzip
import numpy as np
import json


def create_profile(profile_file):
    """
    funkcja do tworzenia slownika ST na podstawie pliku profile sciagnietego z repozytorium enterobase
    :param profile_file: String do pliku profiles, rozne profile (7MLST, cgMLSDT, wgMLST) zawiera rozna liczbe kolumn,
    piwerwsza kolumna to ST, koejne kolumny to informacja jaki wariant w danym miejscu genomu wystepuje (w postaci
     liczby), jesli danego wariantu w danym ST nie ma jest "-".  Pierwszy wiersz to naglowki
    :return: slownik slownikow, jako nadrzedny klucz jest sequence type, kazdemu sequence type odpowiada slownik
     zbudowany z kluczy (nazwy genow z pierwszego naglowka), i jedna wartosc (numerek z pliku)
    """
    slownik_profili = {}
    if not os.path.isfile(profile_file):
        raise IOError('Brak sciezki do pliku z profilem')
    else:
        with open(profile_file) as f:
            for line in f:
                if re.search('ST', line):
                    # jestesmy w wierszy naglowkowym tworzymy liste nazw, ktora posluzy jako klucze w slowniku
                    lista_kluczy = line.rsplit()[1:]  # ST nie jest kluczem
                else:
                    elementy = line.rsplit()
                    # tworzymy slownik dla danego ST
                    slownik_profili[elementy[0]] = {}
                    for indeks, wartosc in enumerate(elementy[1:]):
                        slownik_profili[elementy[0]][lista_kluczy[indeks]] = int(wartosc)
    return slownik_profili, lista_kluczy


def getST(MLSTout, profile_file):
    """
    Funkcja szuka w bazie profili ST o najnizszym dystansie do profilu analizowanej probki. Inpu to profil alleliczny
    probki (MLSTout), oraz lista wszustkich znanych profili w danym schemacie. Odleglosc liczona jest przy pomocy
    funkcji calculate_phiercc_distance, ktora jest reimplementacja funkcji oryginalnie zaimplementowanej w
    programie pHierCC. W przypadku nieobecnosci profilu w znanej bazie identycznego do profilu probki, zwracany
    jest ST o najmniejszej odleglosci, w przypadku wielu takich profili zwracany jest ten ktory ma najnizsze ST.
    :return: lista
    str(ST), min_value, allele_list_lowest_difference, lista_probki, lista_kluczy
    ST - numer najblizszego ST znalezionego w bazie pasujacego do profilu analizowanej probki
    min_value - odleglosc miedzy profilami probk,i a tym znalezionym w bazie
    allele_string_lowest_difference - profil alleliczny dla zwracanego ST, w postaci tab-separarted string
    lista_probki_string - profil alleliczny analizowanej probki w postaci tab-separated string
    lista_kluczy_string - lista zawierajace nazwy loci w danym schemacie w postaci tab separated string
    """

    # pierwszy przelot potrzebujemy listy kluczy/alleli w danym profilu
    with open(profile_file) as f:
        for line in f:
            if re.search('ST', line):
                # jestesmy w wierszy naglowkowym tworzymy liste nazw, ktora posluzy jako klucze w slowniku
                lista_kluczy = line.rsplit()[1:]  # ST nie jest kluczem
                break

    slownik_roznic = {}  # slownik ktory jako klucze ma nazwe profilu, jako wartosc ilosc roznic do "naszej" probki

    # Tworzymy profill alleliczny probki w kolejnosci loci identycznej jak w pliku z baza alleli
    # W przypadku braku danego locus-a dajemy -1
    lista_probki = []
    for klucz in lista_kluczy:
        try:
            lista_probki.append(int(MLSTout[klucz]))
        except KeyError:
            lista_probki.append(-1)

    # inicjalizacja zmiennych
    minimalna_roznica = 3200
    slownik_roznic['local_0'] = minimalna_roznica  # inicjujemy slownik dummy wartoscia
    allele_string_lowest_difference = ''
    lista_probki_string = "\t".join(list(map(str, lista_probki)))
    lista_kluczy_string = "\t".join(list(map(str, lista_kluczy)))

    # iterujemy po pliku z profilami, tym razem liczymy odleglosc miedzy znanymi profilami allelicznymi
    # a profilem probki, przerywamy iteracje po pliku przy znalezieniu profilu z odlegloscia 0
    with open(profile_file) as f:
        for line in f:
            if re.search('ST', line):
                # omijami pierwsza linijke
                pass
            else:
                elementy = line.rsplit()
                # wyciagamy ST  z profilu i allele
                ST_dict = elementy[0]  # ST w pliku z profilami
                lista_ST = list(map(int, elementy[1:]))  # profil alleliczny
                slownik_roznic[ST_dict] = 3200

                slownik_roznic[ST_dict] = calculate_phiercc_distance(lista_probki, lista_ST)

                # updatuje zmienne minimalna roznica i allele_string_lowest_difference
                if slownik_roznic[ST_dict] < minimalna_roznica:
                    minimalna_roznica = slownik_roznic[ST_dict]
                    allele_string_lowest_difference = "\t".join(map(str, elementy[1:]))
                # Przerwyam obliczenia jesli znalazlem odleglosc 0
                if slownik_roznic[ST_dict] == 0:
                    break

    try:
        ST = np.min([int(x) for x, y in slownik_roznic.items() if y == minimalna_roznica])
    except ValueError:
        # obejscie jesli analizujemy baze local gdzie  ST maja postac f'local_{numer}'
        ST = np.min([int(x.split('_')[1]) for x, y in slownik_roznic.items() if y == minimalna_roznica])
        ST = f'local_{ST}'
    return str(ST), minimalna_roznica, allele_string_lowest_difference, lista_probki_string, lista_kluczy_string


def write_novel_sample(profile, output_file):
    """
    Write information about novel sample into a file
    :param profile: String to be saved
    :param output_file: Path to a file where string will be saved, must exists
    :return: bool
    """
    if os.path.exists(output_file):
        with open(output_file, 'a') as f:
            f.write(f'{profile}')
    else:
        raise f'Provided file {output_file} does not exist'


def calculate_phiercc_distance(wektor_x, wektor_y, allowed_missing=0.05):
    """
    Funkcja wyciagnieta z kodu do pHierCC liczy dystans miedzy dwoma profilami
    :param wektor_x: Wektor pierwszego profilu
    :param wektor_y:  Wektor drugiego profilu
    :param allowed_missing: paramter z phierCC uzywany do liczenia dystancu
    :return: int, odleglosc miedzy profilami
    """

    wektor_x = np.array(wektor_x)
    wektor_y = np.array(wektor_y)
    n_loci = len(wektor_y)
    if len(wektor_y) != len(wektor_y):
        raise 'Provided profiles have different lengths'

    rl, ad, al = 0., 1e-4, 1e-4

    # rl ile mamy nie zerowych alleli w drugim analizowanym ST
    ql = np.sum(wektor_x > 0)
    rl = np.sum(wektor_y > 0)
    # ile mamy WSPOLNYCH niezerowych alleli w obu ST
    common_non_zero = (wektor_x > 0) & (wektor_y > 0)
    al += np.sum(common_non_zero)
    # ile mamy NIEZEROWYCH alleli ktore przyjmuja ROZNA wartosc
    ad += np.sum(wektor_x[common_non_zero] != wektor_y[common_non_zero])
    # Wieksza z wartosci ql i rl pomniejszona o dopuszczalna liczbe zerowych alleli
    ll = max(ql, rl) - allowed_missing * n_loci
    ll2 = ql - allowed_missing * n_loci

    if ll2 > al:
        ad += ll2 - al
        al = ll2

    if ll > al:
        ad += ll - al
        al = ll

    dist = int(ad / al * n_loci + 0.5)

    return dist



def parse_MLST_tsv(file_path, long=True, sep="\t"):
    """
    Funkcja do tworzenia slownika MLST na podstawie dwulinijkowego pliku, w ktorym naglowek ma klucze (odpowiadajace
    allelom) a drugiw wiersz to odzzielony ("\t", ";", " ", ",") identyfikator allelu
    :param file_path: string, sciezka do pliku z wynikami z enterobase, chewbacca
    :param long: bool, czy plik zawiera jako pierwsze 2 wpisy identyfikator sekwencji i Sequence type (True), czy tylko
    identyfikator (False)
    :param sep: string, separator kolumn w pliku poanych jako argument file_path
    :return: Funkcja zwraca slownik gdzie luczami sa nazwy alleli a wartosciami wariant
    """
    if long:
        zakres = 2
    else:
        zakres = 1
    slownik_alleli = {}
    with open(file_path, 'r') as f:
        for line in f:
            if re.search('ST', line):
                klucze = line.rsplit(sep)[zakres:]
                klucze = [klucz.split('.')[0] for klucz in klucze]
            else:
                elementy = line.rsplit(sep)[zakres:]
                # tworzymy slownik dla danego ST
                for indeks, wartosc in enumerate(elementy):
                    slownik_alleli[klucze[indeks]] = wartosc

    return slownik_alleli


def parse_MLST_fasta(file_path):
    """
    Funkcja do parsowania wynikow etoki zwracanych jako fasta
    :param file_path: string, sciezka do pliku z wynikami etoki, plik w formacie fasta
    :return: slownik, kluczami sa allele a wartosciami ich wariant
    """
    slownik_alleli = {}
    for record in SeqIO.parse(file_path, "fasta"):
        # opis z etoki jest dosc wystandaryzowany wiec zakladam ze pole description[0] zawiera nazwe wariantu
        # pole 2 zawiera id = wartosc
        # pole 6 zwiera identycznosc sekwencyjna miedzy referencja a tym co jest obserwowane w probce
        try:
            slownik_alleli[record.description.split(' ')[0]] = int(record.description.split(' ')[2].split('=')[1])
        except ValueError:
            # nie znaleziono allelu i funkcja int zwraca blad
            slownik_alleli[record.description.split(' ')[0]] = -1

    return slownik_alleli


def _parse_MLST_blastn(plik):
    """
    Funkcja do parsowania outputu mojego skryptu blastowego w ktorym sa 4 kolumny, pierwsza to nazwa allelu, druga numer
    allelu, 3 to identycznosc sekwenyjna miedzy tym allelem a sekwencja w genomie, 4 to info czy bylo wiele mapowan
    :param plik: str, sciezka do pliku z wynikami funkcji run_blastn_verX.sh
    :return: slownik, kluczami sa nazwy alleli w postaci stringow, warosciami numery alleli w postaci int
    """
    slownik_profili = {}
    if not os.path.isfile(plik):
        raise IOError('Brak sciezki do pliku z profilem')
    with open(plik, 'r') as f:
        for line in f:
            line = line.rsplit()
            try:
                int(line[2])
                if int(line[2]) == 100:
                    slownik_profili[line[0]] = int(line[1])
                else:
                    slownik_profili[line[0]] = -1
            except ValueError:
                slownik_profili[line[0]] = -1
    return slownik_profili


def parse_MLST_blastn(plik):
    """
    Funkcja do parsowania outputu mojego skryptu blastowego sa w nim 2 kolumny klucz i numer allelu (moze byc -1 jesli
    danego allelu nie ma)
    :param plik: str, sciezka do pliku z wynikami funkcji run_blastn_verX.sh
    :return: slownik, kluczami sa nazwy alleli w postaci stringow, warosciami numery alleli w postaci int
    """
    slownik_profili = {}
    if not os.path.isfile(plik):
        raise IOError('Brak sciezki do pliku z profilem')
    with open(plik, 'r') as f:
        for line in f:
            line = line.rsplit()
            try:
                int(line[1])
                slownik_profili[line[0]] = int(line[1])
            except ValueError:
                slownik_profili[line[0]] = -1
    return slownik_profili


def compare_2_allel_dict(slownik1, slownik2, file_prefix='out'):
    """
    Funkcja do porownywania dwoch slownikow alleli w celu okreslenia ich zgodnosci. Funkcja zwraca ilosc kluczy, ktore
    w obu slownikach przyjmuja takie same wartosci. Ponadtow dw a pliki z unikalnymi kluczami dla obu
    slownikow rowniez zapisywane sa do pliko {file_prefix}_slownik1_uniquekeys.txt
    i {file_prefix}_slownik2_uniquekeys.txt

    :param slownik1: slownik
    :param slownik2: slownik, o kluczach identycznych z tymi w slowniku 1
    :param file_prefix: string, funkcja zwraca dwa pliki. Kazdy z nich ma 3 kolumny (nazwa alellu, wersja w slowniku1,
    wersja w slowniku 2). Jeden plik dla kluczy o zgodnych wartosciach miedzy slownikami i dtugi dla kluczy niezgodnych.
    Pliki majaca nazwe {file_prefix}_commoncalue_keys.txt, {file_prefix)_differentvalue_keys.txt
    :return: int, liczba kluczy o zgodnych wartosciach miwedzy slownikami, ponadtwo tworzone sa 2 pliki
    """
    # zrobimy ta petla
    common_keys = 0
    with open(f"{file_prefix}_commoncalue_keys.txt", 'w') as f1, open(f"{file_prefix}_differentvalue_keys.txt", 'w') as f2:
        for allel in slownik1:
            if allel in slownik2.keys():
                if slownik1[allel] == slownik2[allel]:
                    common_keys += 1
                    f1.write(f"{allel}\t{slownik1[allel]}\t{slownik2[allel]}\n")
                else:
                    f2.write(f"{allel}\t{slownik1[allel]}\t{slownik2[allel]}\n")
            else:
                # print(f"Allel {allel} nie wystepuje w drugim slowniku !\n")
                pass
    # na koniec wyrzumy jeszcze do plikow unikalne klucze dla slownikow 1 i 2
    slownik1_keys = set(slownik1.keys())
    slownik2_keys = set(slownik2.keys())
    with open(f'{file_prefix}_slownik1_unique_keys.txt', 'w') as f1, \
            open(f'{file_prefix}_slownik2_unique_keys.txt', 'w') as f2:
        f1.write("\n".join(slownik1_keys - slownik2_keys))
        f2.write("\n".join(slownik2_keys - slownik1_keys))

    return common_keys


def sort_profile(profile_file, sorted_keys_file):
    """
    Funkcja do sortowania/uzupelaniania profili o brakujace wpisy, zgodnie z zewnetrzna lista
    :param profile_file: str, sciezka do pliku z profile. Plik ma dwa wiersze, perwszy to lista loci (uwaga niezaleznie
    od nazwy pierwsza kolumna ZAWSZE traktowana jest jako nazwa dla Sequence type), drugi wiersz to int-y z allelami
    :param sorted_keys_file: sciezka do pliku z jednym wierszem gdzie "\t" oddzielono allele. Uwaga NIE podawac kolumny
    identyfikujacej ST
    :return: plik o nazwie basename(profile_file)_sorted_allells.txt
    """

    slownik_alleli = {}
    out_name = f"{os.path.basename(profile_file).split('.')[0]}_sorted_allells.txt"

    sorted_keys = open(sorted_keys_file).readlines()
    sorted_keys = sorted_keys[0].rsplit()

    with open(profile_file) as f:
        for line in f:
            if re.search('ST', line):
                klucze = line.rsplit()
                # czasami jako output jest podawany nazwa.fasta (od nazwy pliku do szukania alleli)
                klucze = [klucz.split('.')[0] for klucz in klucze]
            else:
                elementy = line.rsplit()
                for indeks, wartosc in enumerate(elementy):
                    slownik_alleli[klucze[indeks]] = wartosc

    with open(out_name, 'w') as f:
        out_msg = '800000'  # dummy identyfikator ?
        for klucz in sorted_keys:
            try:
                slownik_alleli[klucz]
            except KeyError:
                slownik_alleli[klucz] = '-1'
            out_msg = out_msg + '\t' + slownik_alleli[klucz]

        f.write("ST" + "\t" + "\t".join(sorted_keys) + "\n")
        f.write(out_msg + "\n")

    return True


if __name__ == '__main__':
    pass
