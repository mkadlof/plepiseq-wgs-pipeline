#!/usr/bin/env python3

'''
prosty konwerter outputu varscan do vcf-a. jedyny argument to plik stowrzony przez varscan
dokumantacja outputu varscan jest tu http://varscan.sourceforge.net/using-varscan.html
UPDATE 16.03
w toku EQA zauwazylem ze moj skrypt nie radzi sobie z sytuacja w varscan 
MN908947.3      22193   A       */-ATT  195     403
gdzie powinno byc
MN908947.3      22193   .       AATT    A
a jest 
ESIB_EQA_2023.SARS2.16/output_detected_variants_varscan_oryg.vcf:MN908947.3     22193   .       *       A

nie wiem czemu czasami varscan daje wniosek jakby nasza probka byla diploidem i mogla byc heterozygotyczna, ale poprzednie casey zawsze
dawaly sytuacje ze heterozygota miala  akie same wersje alleelu
tak wiec wprowadzamy malo poprawke do skryptu
jesli jest +/- i sample jest hetwero to:
    rozdzielamy alt allele na dwa split "/"
    potem porownujemy czy sa identyczne 
    jak nie to czy ktorys jest gwiazdka/kropka jaeli tak to wypieramy ten drugi
    jesli sa rozne to wybieramy dluzsza insercje/delecje
    printujemy komunikat ze takie sytuacje zaszly na out-a


update 22 06
varscan ma dosc liberalne podejscie do nukleotydow ambigous nawet 70% do 30% stosunek to daje, robimy wiec tak,
ze  ALT ktory jest w pozycji 19 w wierszu zastapi ten ambigous kolumna 4 jesli jego pokrycie jest powyzej 55 % kolumna 7 (uwaga ma procenty wiec trzeba to przerobic
na floata
'''

import sys
import os
from datetime import datetime

nazwa_pliku = os.path.splitext(os.path.basename(sys.argv[1]))[0]
current_date = datetime.now().strftime("%Y%m%d")

upper_ambig=float(sys.argv[2])
pval = float(sys.argv[3])

with open(f'{nazwa_pliku}.vcf', 'w') as f:
    #dodajmy naglowek vcf-a
    f.write('##fileformat=VCFv4.0\n')
    f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    f.write(f'##fileDate={current_date}\n')
    f.write('##source=varscan\n')
    f.write('##reference=User-provided\n')
    f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">\n') #Read1 + Read2
    f.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n') #  VarFreq
    f.write('##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse,'
            ' alt-forward and alt-reverse bases">\n') # Reads1Plus,  Reads1Minus, Reads2Plus, Reads2Minus
    f.write('##INFO=<ID=QR,Number=1,Type=Integer,Description="Average map quality of ref reads">\n')  # MapQual1
    f.write('##INFO=<ID=QA,Number=1,Type=Integer,Description="Average map quality of alt reads">\n')  # MapQual2
    f.write('##INFO=<ID=GL,Number=1,Type=Float,Description="Significance of variant read count vs. expected baseline erro">\n')  # Pvalue
    f.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    with open(sys.argv[1]) as f1:
        i = 0
        for line in f1:
            if i == 0:
                i += 1
                continue

            line = line.split('\t')
            CHROM=line[0]
            POS=line[1]
            ID ='.'
            if len(line[3]) > 1:
                # INDEL
                if '+' in line[3]:
                    # insercja np
                    # N908947.3 28262 G +AACA/+AACA
                    # musi byc poprawiona na
                    # MN908947.3 28262   .       G       GAACA
                    # czyli w miejscu '+' wstawimy G
                    REF=line[2]
                    if '/' in line[3]:
                        left,right = line[3].split('/')
                        if left != right:
                            print(f'Na pozycji {POS} jest heterozygota {left}/{right}')
                            if left != "*" and len(left) >= len(right):
                                # ten warunek na gwiazdke nie ma sensu bo * raczej nigdy nie bedzie dluzsza
                                ALT=line[3].replace('+', line[2]).split('/')[0]
                            elif right != "*" and len(right) > len(left):
                                ALT=line[3].replace('+', line[2]).split('/')[1]
                        else:
                            # Stara wersja
                            #REF=line[2]
                            ALT=line[3].replace('+', line[2]).split('/')[0]
                            
                    else:
                        # Z alt wywalilem tego splita line[3].replace('-', line[2]).split('/')[0]
                        #REF=line[2]
                        ALT=line[3].replace('+', line[2])
                elif '-' in line [3]:
                    # delecja np.
                    # MN908947.3 11287 G  -TCTGGTTTT/-TCTGGTTTT
                    # musi byc poprawiona na
                    # MN908947.3 11287 . GTCTGGTTTT G
                    # czyli '-' zamieniamy na to co jest w polu REF. Na koncu ALT staje sie REF, a REF ALT
                    ALT = line[2]
                    if '/' in line[3]:
                        pass
                        left,right = line[3].split('/')
                        if left != right:
                            print(f'Na pozycji {POS} jest heterozygota {left}/{right}')
                            if left != "*" and len(left) >= len(right):
                                REF= line[3].replace('-', line[2]).split('/')[0]
                            elif right != "*" and len(right) > len(left):
                                REF= line[3].replace('-', line[2]).split('/')[1]
                        else:
                             REF= line[3].replace('-', line[2]).split('/')[0]
                    else:
                        REF= line[3].replace('-', line[2]).split('/')
                    
                    #ALT = line[2]
                else:
                    raise Exception(f'Unsuppoted scenario either + or - must be present in vcf!\n{line}')
            else:
                #SNP
                REF = line[2]
                if float(line[6].split('%')[0])/100 > upper_ambig :
                    ALT = line[18].split()[0]
                else:
                    ALT = line[3]
            # program nie ma quality jako takiego ale ma p-value 0 - oznacza swietny wariant , dajmy sobie threshols hmmm 0.01 ? jesli taka
            # bedzie wartosc damy quality 30 jak nizej to qual 1
            if float(line[11]) <= pval:
                QUAL = 30 # za quality bedzie p-value
            else:
                QUAL=1
            FILTER = 'PASS'
            INFO=f'DP={int(line[4]) + int(line[5])};AF={float(line[6].strip("%"))/100};DP4={line[14]},{line[15]},{line[16]},{line[17]};' \
                 f'QR={line[9]};QA={line[10]};GL={QUAL}'
            f.write(f'{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\n')



