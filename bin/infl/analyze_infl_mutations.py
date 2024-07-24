#!/usr/bin/env python3

"""
Compared to the version with Nanopore, we have additional variables for
file locations (Output_dir) and data (data_dir).
"""

import os
import re
import shutil
import subprocess
import sys
from tempfile import NamedTemporaryFile
from typing import Dict

import numpy as np
from Bio import SeqIO


# Przykladowe wywolanie w kontenerze z illuminy (gdzie jest skopiowany i poprawion w .opt/docker. repo netclade)
# to 0.1 to threshold na ilsc x w sekwencji
# python analyze_Influenza_mutations.py  ESIB_EQA_2023.INFL4.01.fasta 0.1

# 1 identyfikujemy skad jest nasz wirus na podstawie sekwencji HA

def determine_subtype(subtype):
    """
    :param subtype:
    :return:
    """

    # ta sekwencja lepeij zeby zawierala tylko HA

    # zwracamy nazwy tak by pasowaly do sciezek z /opt/docker/nextclade/data/flu (konetner illumina do sars-a ! poki co'

    if 'H1N1' in subtype and 'pdm09' in subtype:
        return 'h1n1pdm'
    elif 'H1N1' in subtype and 'pdm09' not in subtype:
        return 'h1n1'
    elif 'H3N2' in subtype:
        return 'h3n2'
    elif 'B/Yam' in subtype:
        return 'yam'
    elif 'B/Vic' in subtype:
        return 'vic'
    elif 'unknown' in subtype:
        raise Exception(f'unknown subtype')
    else:
        raise Exception(f'Unknown case {subtype}')


def check_protein_quality(plik):
    """'
    Sprawdzamy obecnosc X w sekwencji. Jesli sumarycznie jest wiecej niz threshold, definioany
    jako ich frakcja sekwencja dostaje flage ze jest zla. Uwaga sekwencja jest alignowana wzgledem refencji
    wiec usuawamy '-' do liczenia tego stosunku
    :param plik:
    :return:
    """
    tekst = 'Brak sekwencji'
    stosunek = 1
    record = SeqIO.parse(plik, 'fasta')
    for r in record:
        ilosc = int(r.seq.count('X'))
        dlugosc = float(len(r.seq.replace('-', '')))
        stosunek = float("%.4f" % (ilosc / dlugosc))
        tekst = f'Sekwencja {r.id} ma dlugosc:\t{dlugosc}\tIlosc X:\t{ilosc} {stosunek * 100}%\n'
        return tekst, stosunek
    return tekst, stosunek


def check_genome_segment_quality(plik):
    """'
    Sprawdzamy obecnosc X w sekwencji. Jesli sumarycznie jest wiecej niz threshold, definioany
    jako ich frakcja sekwencja dostaje flage ze jest zla. Uwaga sekwencja jest alignowana wzgledem refencji
    wiec usuawamy '-' do liczenia tego stosunku
    :param plik:
    :return:
    """
    tekst = 'Brak sekwencji'
    stosunek = 1
    record = SeqIO.parse(plik, 'fasta')
    for r in record:
        ilosc = int(r.seq.count('N'))
        dlugosc = float(len(r.seq.replace('-', '')))
        stosunek = float("%.4f" % (ilosc / dlugosc))
        tekst = f'Sekwencja {r.id} ma dlugosc:\t{dlugosc}\tIlosc X:\t{ilosc} {stosunek * 100}%\n'
        return tekst, stosunek
    return tekst, stosunek


def get_fasta_header(input_fasta):
    record = SeqIO.parse(input_fasta, 'fasta')
    for r in record:
        return f'>{r.description}'


def get_fasta(input_fasta, sample, what):
    """
    :param input_fasta:
    :param sample:
    :param what:
    :return:
    """
    record = SeqIO.parse(input_fasta, 'fasta')
    two_HA = 0
    for r in record:
        if what in r.id:
            if two_HA > 0:
                raise Exception(f'File {input_fasta} contains more than one fragments with {what} in its name, quiting !')
            # cala zabawa z okreslaniem subtypu
            with NamedTemporaryFile(dir='.', prefix=f'{sample}_', suffix=f'_{what}.fasta', mode='w', delete=False) as file:
                file.write(f'>{r.id}\n')
                file.write(f'{str(r.seq)}\n')
            two_HA += 1
    # gdyby nie bylo danego bialka w wynikach z jakiegos powodu tworzymy dummy file
    # ktory nie bedzie analizowany
    if two_HA == 0:
        with NamedTemporaryFile(dir='.', prefix=f'{sample}_', suffix=f'_{what}.fasta', mode='w', delete=False) as file:
            file.write(f'>{what}\n')
            file.write(f'XXX\n')
    return file


def _run_infinity(file, input_fasta):
    """
    Funkcja do puszczania infinity programu okreslajacego clad do ktorego nalezy podana sekwencja HA.
    Aktualnie obsolete zastapione tylko nextclade, choc w przyszlosci moze nada sie do avian ?
    :param file:
    :param input_fasta:
    :return:
    """
    if 'run_infinity.R' not in os.listdir('.'):
        raise Exception('Pamietaj o umieszczeni uskryptu run_infinity.R w katalogu gdzie wywolujesz skrypt!')
    komenda = f'Rscript run_infinity.R {file.name}'
    infinity = subprocess.Popen(komenda, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in infinity.stderr:
        line = line.decode('UTF-8').rstrip()
        if 'Error' in line:
            raise Exception(f"Infinity could not process provided sequence from {input_fasta}")

    for line in infinity.stdout:
        line = line.decode('UTF-8').rstrip()
        if 'Clade' in line:
            continue
        print(line)
        subtype = determine_subtype(line.split(';')[2])
    os.remove(file.name)
    return subtype


# 2 korzystajac z odpowiedniej bazy okreslamy aminokwasowa sekwencje danego HA
def get_protein_sequence(segment_fasta, referene_fasta, gff_file, what, out_dir, qc_threshold, typ, sample_name):
    # referencyjna sekwencja
    # uwaga w /opd.docker/nextclade/data/flu itd .. w genmap musialem recznie poprawic opis z gene_mamae "ALA" na gene_name="ALA" !!!!
    # zrobilem to tylko dla HA, NA i PA bo te geny sa wazne dla opornosci
    if f'Reference_{what}_{typ}.fasta' not in os.listdir('.'):
        komenda = f'nextalign run -r {referene_fasta} -m {gff_file} -g {what} -O {out_dir} -n {what} {referene_fasta}'
        p = subprocess.Popen(komenda, shell=True, stdout=subprocess.PIPE)
        p.communicate()
        shutil.copyfile(f'./{out_dir}/{what}_gene_{what}.translation.fasta', f'Reference_{what}_{typ}.fasta')
        shutil.rmtree(f'{out_dir}')

    # teraz sekwencja z probki
    with open(f'QC_report_{what}.txt', 'w') as f:
        komenda = f'nextalign run -r {referene_fasta} -m {gff_file} -g {what} -O {out_dir} -n {what} {segment_fasta}'
        p = subprocess.Popen(komenda, shell=True, stdout=subprocess.PIPE)
        p.communicate()
        os.remove(f'{segment_fasta}')
        tekst, stosunek = check_protein_quality(f'{out_dir}/{what}_gene_{what}.translation.fasta')
        f.write(tekst)

        if stosunek < qc_threshold:
            shutil.copyfile(f'{out_dir}/{what}_gene_{what}.translation.fasta', f'{sample_name}_{what}.fasta')
        shutil.rmtree(f'{out_dir}')

    return True


# 4 alignment na odpowiednia referencje

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


# 5 mutacje z odpowiedniej tabelki

def determine_muation(ref_seq, target_seq):
    """
    Dwa stringi zawierajace A-Z i "-" lub "N" (ale tylko target), mamy nastepujace stany,
    1. brak mutacji
    2. mutacja roznica charakteru miedzy referencja a targetem
    3. rozpoczecie wstawki
    4. kontynuacja wstawki
    5. rozpoczecie insercji
    6. kontynuacja insercji
    raportujemy jako
    what:R4L punktowa
    R4_del2 - delecja 2 aminokwasow w target po R4
    R4_ins4 - insercja 4 aminokwasow w target po R4
    :param ref_seq:
    :param target_seq:
    :return:
    """

    # sekwencje sa identycznej dlugosci wiec

    numeracje_ref = 1  # od 1
    previous_ref = 'X'
    previous_target = "X"
    poczatek_delecji = 0
    dlugosc_insercji = 0
    poczatek_insercji = 0
    dlugosc_delecji = 0

    lista_mutacji = []

    for i in range(len(ref_seq)):

        # sprawdz czy nie skonczyles insercji
        if dlugosc_insercji > 0 and ref_seq[i] != '-':
            lista_mutacji.append(f'{ref_seq[poczatek_insercji - 1]}{poczatek_insercji}ins{dlugosc_insercji}')
            dlugosc_insercji = 0
            poczatek_insercji = 0
        # sprawdz czy nie skonczyles delecji

        if dlugosc_delecji > 0 and target_seq[i] != '-':
            # lista_mutacji.append(f'{ref_seq[poczatek_delecji - 1]}{poczatek_delecji}del{dlugosc_delecji}')
            lista_mutacji.append(f'del_{poczatek_delecji}_{poczatek_delecji + dlugosc_delecji}')
            dlugosc_delecji = 0
            poczatek_delecji = 0

        if ref_seq[i] == target_seq[i]:
            pass
        elif ref_seq[i] != target_seq[i] and ref_seq[i] != '-' and target_seq[i] != '-':
            # prawdiwa mutacja
            lista_mutacji.append(f'{ref_seq[i]}{numeracje_ref}{target_seq[i]}')
        elif ref_seq[i] != target_seq[i] and ref_seq[i] == '-' and previous_ref != '-':
            # zaczynamy insercje
            poczatek_insercji = numeracje_ref
            dlugosc_insercji = 1
        elif ref_seq[i] != target_seq[i] and ref_seq[i] == '-' and previous_ref == '-':
            dlugosc_insercji += 1
        elif ref_seq[i] != target_seq[i] and target_seq[i] == '-' and previous_target != '-':
            # poczatek delecji
            poczatek_delecji = numeracje_ref
            dlugosc_delecji = 1
        elif ref_seq[i] != target_seq[i] and target_seq[i] == '-' and previous_target == '-':
            dlugosc_delecji += 1
        else:
            print(f'Unknown case na pozycji alignmentu {i} ref {ref_seq[(i - 3):(i + 1)]}target {target_seq[(i - 3):(i + 1)]} ')

        if ref_seq[i] != '-':
            numeracje_ref += 1

        previous_ref = ref_seq[i]
        previous_target = target_seq[i]

    if dlugosc_insercji > 0 and ref_seq[i] != '-':
        lista_mutacji.append(f'{ref_seq[poczatek_insercji - 1]}{poczatek_insercji}ins{dlugosc_insercji}')
        dlugosc_insercji = 0
        poczatek_insercji = 0
    # sprawdz czy nie skonczyles delecji

    if dlugosc_delecji > 0 and target_seq[i] != '-':
        lista_mutacji.append(f'{ref_seq[poczatek_delecji - 1]}{poczatek_delecji}ins{dlugosc_delecji}')
        dlugosc_delecji = 0
        poczatek_delecji = 0

    return lista_mutacji


def higest_rest(tekst):
    if 'HRI' in tekst:
        return 'HRI'
    elif 'RI' in tekst:
        return 'RI'
    elif 'NI' in tekst:
        return 'NI'
    elif '?' in tekst:
        # nie badano, zakladamy normal inhibition
        return 'NI'
    else:
        return 'Err'


def parse_resistance_NA(plik):
    """
    Parser do pliku z opornosci jest 5 kolumn pierwsza to mutacja (badz kombinacja mutacji oddzielonych +, lub delecja_od_do)
    kolejne 4 kolumy to tekstowy opis opornosci kolejnych 4 lekow:
    Oseltamivir, Zanamivir, Peramivir, Laninamivir
    Na razie wystarczy wiedziec czy jest tam NI, RI, HRI
    :param plik:
    :return:
    """
    slownik_opornosci_oseltamivir = {}
    slownik_opornosci_zanamivir = {}
    slownik_opornosci_peramivir = {}
    slownik_opornosci_laninamivir = {}
    with open(plik) as f:
        for line in f:
            if 'Mutation' in line:
                continue
            try:
                klucz, oselt, zana, pera, lana = line.split('\t')
                slownik_opornosci_oseltamivir[klucz] = higest_rest(oselt)
                slownik_opornosci_zanamivir[klucz] = higest_rest(zana)
                slownik_opornosci_peramivir[klucz] = higest_rest(pera)
                slownik_opornosci_laninamivir[klucz] = higest_rest(lana)
            except:
                print(f'Error for {line}')
    return slownik_opornosci_oseltamivir, slownik_opornosci_zanamivir, slownik_opornosci_peramivir, slownik_opornosci_laninamivir


def parse_resistance_PA(plik):
    """
    Parser do pliku z opornosci, sa 2 kolumny pierwsza to mutacja (tylko punktowe) opornosc na baloxavir jest kolumnie 2
    Na razie wystarczy wiedziec czy jest tam NI, RI, HRI. W odroznieniu od NA tu NIE ma naglowka
    :param plik:
    :return:
    """
    slownik_opornosci = {}

    with open(plik) as f:
        for line in f:
            if 'Mutation' in line:
                continue
            try:
                klucz, balo = line.split('\t')
                slownik_opornosci[klucz] = higest_rest(balo)
            except:
                print(f'Error for {line}')

    return slownik_opornosci


def get_sample_status(lista_mutacji, slownik_opornosci):
    mutacje_opornosciowe = {}
    for klucz in slownik_opornosci:
        status = slownik_opornosci[klucz]
        if status in ['RI', 'HRI']:
            if set(klucz.split('+')).issubset(set(lista_mutacji)):
                mutacje_opornosciowe[klucz] = status
        else:
            # mutacje nieistotne dla danego leku, nie raportuje raczej bo po co
            continue

    status = 'S'
    # teraz okreslamy jaka jest najmocniejsza mutacja i zwracamy stringa "S" dla NI, I dla RI i R dla HRI i 'U' jesli nie ma opornosci
    for klucz in mutacje_opornosciowe:
        if mutacje_opornosciowe[klucz] == 'RI' and status in ['S']:
            status = 'I'
        elif mutacje_opornosciowe[klucz] == 'HRI' and status in ['S', 'I']:
            status = 'R'

    return mutacje_opornosciowe, status


if __name__ == '__main__':
    # Modyfikacja skryptu z 2023
    # inputem jest fasta z bialkami (NA i PA)
    # to jest wyolanie glownej funkcji dla pojedynczego sample'a
    plik_fasta = sys.argv[1]  # plik z sekwencjami bialek, moze byc multifasta
    podtyp = sys.argv[2]  # podtyp analizowanego wirusa (obslugujem H1N1, H3N2, H5N1, H7N9, Yamagata, Victoria)
    sample_name = sys.argv[3]  # nazwa sample uzywana bedzie w podplikach tworzonych wprzez skrypt
    data_dir = sys.argv[4]  # sciezka gdzie trzymane sa dane w kontenerze
    output_dir = sys.argv[5]  # sciezka w ktorym katalogu zapisywac wyniki

    # Definicja zmiennych
    lista_mutacji_NA_opornosc_oselt = ['UNK']
    lista_mutacji_NA_opornosc_zana = ['UNK']
    lista_mutacji_PA_opornosc_balo = ['UNK']
    status_oselt = 'UNK'
    status_zana = 'UNK'
    status_balo = 'UNK'

    # plik* to obiekty clasy templfile ich nazwa jest w obiek.name, trzymaja genomiczna sekwencje danego segmentu
    plik_referencja_NA = f"/{data_dir}/resistance/{podtyp}/{podtyp}_NA_reference.fasta"
    opornosci_NA = f"{data_dir}/resistance/{podtyp}/{podtyp}_NA_resistance.txt"
    plik_fasta_NA = get_fasta(input_fasta=plik_fasta, what='NA', sample=sample_name)
    tekst_NA_genome, stosunek_NA_genome = check_protein_quality(plik_fasta_NA.name)
    if stosunek_NA_genome > 0.2:
        print('Sekwencja bialka NA zawiera ponad 20% N-ek, nie analizuje')

    else:
        # alignment
        NA_alignment = align_fasta(plik_fasta_NA.name, plik_referencja_NA)
        NA_ref_header = get_fasta_header(plik_referencja_NA)
        NA_target_header = get_fasta_header(plik_fasta_NA.name)
        if NA_ref_header == NA_target_header:
            raise Exception(f'Sekwencje NA maja taki sam header, wychodze')
        # To pisze w miejscu wywolania i po pierwsze tworzy liste mutacji
        # w postaci {aminokwas_referencyjny}_{numer w referncji}_{aminokwas w target}
        # delecje i insercje opisane sa w postaci f'del_{poczatek_delecji}_{poczatek_delecji+dlugosc_delecji}'
        # i f'{ref_seq[poczatek_insercji - 1]}{poczatek_insercji}ins{dlugosc_insercji}'

        with open(f'{output_dir}/Mutation_report_NA.txt', 'w') as f:
            lista_mutacji_NA = determine_muation(NA_alignment[NA_ref_header], NA_alignment[NA_target_header])
            f.write(f'{sample_name}\t{";".join(lista_mutacji_NA)}\n')
        # slowniki opornosci
        slownik_opornosci_oseltamivir_NA, \
            slownik_opornosci_zanamivir_NA, \
            slownik_opornosci_peramivir_NA, \
            slownik_opornosci_laninamivir_NA = parse_resistance_NA(opornosci_NA)

        # wyniki
        lista_mutacji_NA_opornosc_oselt, status_oselt = get_sample_status(lista_mutacji_NA, slownik_opornosci_oseltamivir_NA)
        lista_mutacji_NA_opornosc_zana, status_zana = get_sample_status(lista_mutacji_NA, slownik_opornosci_zanamivir_NA)
        lista_mutacji_PA_opornosc_balo, status_balo = ['UNK'], 'UNK'

    if podtyp in ['H1N1', 'H3N2', 'Yamagata', 'Victoria']:
        # opornosci na PA sa znane tylko dala tych podtypow
        plik_referencja_PA = f"/{data_dir}/resistance/{podtyp}/{podtyp}_PA_reference.fasta"
        opornosci_PA = f"/{data_dir}/resistance/{podtyp}/{podtyp}_PA_resistance.txt"
        plik_fasta_PA = get_fasta(input_fasta=plik_fasta, what='PA', sample=sample_name)
        tekst_PA_genome, stosunek_PA_genome = check_protein_quality(plik_fasta_PA.name)
        if stosunek_PA_genome > 0.2:
            print('Bialko PA zawiera ponad 20% N-ek, nie analizuje')
        else:
            PA_alignment = align_fasta(plik_fasta_PA.name, plik_referencja_PA)
            PA_ref_header = get_fasta_header(plik_referencja_PA)
            PA_target_header = get_fasta_header(plik_fasta_PA.name)
            if PA_ref_header == PA_target_header:
                raise Exception(f'Sekwencje PA maja taki sam header, wychodze')

            with open(f'{output_dir}/Mutation_report_PA.txt', 'w') as f:
                lista_mutacji_PA = determine_muation(PA_alignment[PA_ref_header], PA_alignment[PA_target_header])
                f.write(f'{sample_name}\t{";".join(lista_mutacji_PA)}\n')

            slownik_opornosci_baloxavir_PA = parse_resistance_PA(opornosci_PA)
            lista_mutacji_PA_opornosc_balo, status_balo = get_sample_status(lista_mutacji_PA,
                                                                            slownik_opornosci_baloxavir_PA)

    # zapisywanie raportow
    with open(f'{output_dir}/Mutation_report_only_resistance.txt', 'w') as f:
        f.write(f'{sample_name}\t{";".join(lista_mutacji_NA_opornosc_oselt)}\t{";".join(lista_mutacji_NA_opornosc_zana)}\t{";".join(lista_mutacji_PA_opornosc_balo)}\n')

    with open(f'{output_dir}/Output_for_EQA.txt', 'w') as f:
        f.write(f'{plik_fasta}\t{status_oselt}\t{status_zana}\t{status_balo}\n')
