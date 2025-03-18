#!/usr/bin/env python3
'''
Skrypt zamienia SNP-y wykryte przez medaka na ambigous jesli wedlug varscana uzycie allelu alternatynwego jest ponizej upper_ambig
Nie testujemy dolnej granicy, bo medaka ma zawsze racje i lower ambi traci tu znaczenie
wersja na podstawie _5, usunalem tylko stary/nieuzywany kod
'''

import sys
import vcf

artic_vcf = sys.argv[1]
varscan_vcf = sys.argv[2]

upper_ambig=float(sys.argv[3])
out = sys.argv[4]
spare = sys.argv[5] # artic vcf de facto medaka po annotacji i filtrowaniu moze byc pusta w drugiej rundzie w takim wypadu
# nie mam pola record z vcf-a i program nie wprowadza pozycji ambigous 
# dlatego biore nie pofiltroany plik ambigous, gorzej co jesli i ten jest pusty ... 

cov = int (sys.argv[6])

lines_from_varscan = {}
with open(varscan_vcf) as f:
    for line in f:
        if 'Chrom' in line:
            pass
        else:
            line = line.split('\t')
            CHROM = line[0]
            if CHROM not in lines_from_varscan.keys():
                lines_from_varscan[CHROM] = {}
            POS = int(line[1])
            COVERAGE=int(line[14]) + int(line[15]) + int(line[16]) + int(line[17])
            if len(line[3]) > 1:
                # INDEL-e
                if '+' in line[3]:
                    pass
                elif '-' in line[3]:
                    pass
                else:
                    pass
            else:
                # SNP-y
                REF = line[2]
                if float(line[6].split('%')[0])/100 > upper_ambig:
                    pass
                else:
                    if COVERAGE > cov:
                        ALT = line[3]
                        lines_from_varscan[CHROM][POS] = line


# na tym etapie w slowniku mamy mutacje ktore sa wedlug varscana ambigous
# jako wartosc slownik ma po prostu linijke z pliku txt

artic_vcf_reader = vcf.Reader(filename=artic_vcf)
vcf_writer = vcf.Writer(open(out, 'w'), artic_vcf_reader)

# do ostatecznego pliku zapisujemy tylko te pozycje ktore sa w medace
# jedyne na co pozwalamy to dla pozycji ktore sa w medace i w slowniku z varscana
# podmienic symbol dla allelu alternatywngo

for record in artic_vcf_reader:
    try:
        # w drugiej rundzie czesc segmentow grypy moze nie miec mutacji i moze ich brakowac
        lines_from_varscan[record.CHROM]
    except:
        lines_from_varscan[record.CHROM] = {}

    if record.POS not in lines_from_varscan[record.CHROM].keys():
        print(f'{record.CHROM}\t{record.POS} is not ambigous')
        vcf_writer.write_record(record)
    else:
        if record.var_type == 'snp':
            wiersz=lines_from_varscan[record.CHROM][record.POS]
            print(f'{record.CHROM}\t{record.POS} is ambigous')
            QUAL = 30 # varscan nie ma walidowac medaki dajemy qual 30
            # tworzymy record tak jak 3ma go pyvcf3
            # pola INGO, format, samples i sample_index kopiujemy z vcf-a medaki tak by zgadzal sie ich format
            a = vcf.model._Record(CHROM = wiersz[0], POS = int(wiersz[1]), ID = '',  REF = wiersz[2], ALT = [vcf.model._Substitution(nucleotides = wiersz[3])], QUAL = QUAL,\
                FILTER = [], INFO = record.INFO, FORMAT = record.FORMAT, samples = record.samples, sample_indexes = record._sample_indexes)
            vcf_writer.write_record(a)
            del(lines_from_varscan[record.CHROM][record.POS])
        else:
            print(f'{record.CHROM}\t{record.POS} is ambigous but not snp')
            vcf_writer.write_record(record)
            del(lines_from_varscan[record.CHROM][record.POS])

try:
    record.samples
except:
    # nie ma record-u zadnego bierzemy pole record z pliku annotated przed filtrowaniem
    # ale zmieniamy w nim pola GQ i DP, te wartosci przyjma wartosc 100 
    artic_vcf_reader = vcf.Reader(filename=spare)
    for record in artic_vcf_reader:
        CallDataTuple = vcf.model.make_calldata_tuple(['GT', 'GQ'])
        record.samples=[vcf.model._Call(sample='SAMPLE', site=record.genotype('SAMPLE').site, data=CallDataTuple('1','100'))]
        record.INFO['DP'] = 100
        break


print('Medaka nie znalazla nastepujacych pozycji ambigous wlaczam je do analizy jesli maja poprawne p-value')
print(lines_from_varscan)
for CHROM in lines_from_varscan:
    for pozycja, wiersz in lines_from_varscan[CHROM].items():
        QUAL = 30 # za quality bedzie p-value
        a = vcf.model._Record(CHROM = wiersz[0], \
                POS = int(wiersz[1]), \
                ID = '.', \
                REF = wiersz[2], \
                ALT = [vcf.model._Substitution(nucleotides = wiersz[3])], # ALT jest lista klasy ._Substitution \
                QUAL = QUAL, \
                FILTER = [], \
                INFO = record.INFO, FORMAT = record.FORMAT, samples = record.samples, sample_indexes = record._sample_indexes)
        vcf_writer.write_record(a)
