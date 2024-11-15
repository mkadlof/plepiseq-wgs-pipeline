#!/usr/bin/env python3
'''
Ten skrypt dodaje do medaki wybrane wyniki varscana,  w tym wypadku te SNP-y ! ktorych uzycie allelu alternatywnego jest wiekjsze 
niz 80% a wyniki nie znalazly sie w z jakis powodow w medace (patrze dokumentacja i opisy EQA). Filtrowanie odbywa sie prostym awkiem 
ten program po prostu formatuje plik .txt varscana do postaci vcf-a jak medaka
'''

import sys
import vcf

artic_vcf = sys.argv[1] # output medaki .vcf.gz
varscan_vcf = sys.argv[2] # output varscana .txt
out = sys.argv[3] # output programu .vcf

# nie mam pola record z vcf-a i program nie wprowadza pozycji ambigous 
# dlatego biore nie pofiltroany plik ambigous, gorzej co jesli i ten jest pusty ... 

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
        print(f'{record.CHROM}\t{record.POS} is not present in varscan')
        vcf_writer.write_record(record)
    else:
        wiersz=lines_from_varscan[record.CHROM][record.POS]
        print(f'{record.CHROM}\t{record.POS} is present in varscan')
        vcf_writer.write_record(record)
        del(lines_from_varscan[record.CHROM][record.POS])


# Tworzymy pusty record ktorego bedziemy uzywali aby moc skonwertowac wyniki z pliku .txt varscana na format .vcf z medaki /wartosci jak GQ, DP itd beda mialy wartosc 100/
artic_vcf_reader = vcf.Reader(filename=artic_vcf)
for record in artic_vcf_reader:
    CallDataTuple = vcf.model.make_calldata_tuple(['GT', 'GQ'])
    record.samples=[vcf.model._Call(sample='SAMPLE', site=record.genotype('SAMPLE').site, data=CallDataTuple('1','100'))]
    record.INFO['DP'] = 100
    record.INFO['DPS'] = '100,100'
    break


print('Medaka nie znalazla nastepujacych pozycji laczam je do analizy')
print(lines_from_varscan)
for CHROM in lines_from_varscan:
    for pozycja, wiersz in lines_from_varscan[CHROM].items():
        QUAL = 100 # za quality bedzie p-value
        a = vcf.model._Record(CHROM = wiersz[0], \
                POS = int(wiersz[1]), \
                ID = '.', \
                REF = wiersz[2], \
                ALT = [vcf.model._Substitution(nucleotides = wiersz[3])], # ALT jest lista klasy ._Substitution \
                QUAL = QUAL, \
                FILTER = [], \
                INFO = record.INFO, FORMAT = record.FORMAT, samples = record.samples, sample_indexes = record._sample_indexes)
        vcf_writer.write_record(a)
