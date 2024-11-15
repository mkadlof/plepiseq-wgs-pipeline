#!/usr/bin/env python3
"""
This script is used to get coverage at a given position(s)
and usage of alternative allele(s)

for most simple cases
chrom 10 A T it is relatively straightforward to calculate coverage and usage of T
for cases
chrom 10 TAT AAA we caclulate coverage at EACH position and usage at EACH position and take a mean
for deletions and insertions
chrom 16 C CTTTA we take coverage and position 16 and calculate how many reads have is_ins characteristi
same for deletion (is_del characteristic)

"""
import sys
import pysam
import vcf
import numpy as np

# okreslamy czy mamy do czynienia z delecja, insercja, "dlugim" snp-em czy pojedycznym SNPem


def deterine_status(segment_name, position_start, ref_allel, alt_allel):
    if len(ref_allel) == len(alt_allel) and len(ref_allel) == 1:
        status = "snp"
    elif len(ref_allel) > len(alt_allel):
        status = "del"
    elif len(ref_allel) < len(alt_allel):
        status = "ins"
    else:
        raise Exception(f'Unknown case for {segment_name}\t{position_start}\t{ref_allel}\t{alt_allel}')
    return status


def read_vcf(vcf_file):
    lines_from_vcf = {}  # slownik nadrzedny klucz to contig, podrzedny to pozycja np lines_from_vcf['chr1'][134]
    # wartoscia bedzie lista dwuelementow ['AAA', 'A', "short"]  gdzie pierwszy element to wartosc REF a drugi ALT, a
    # trzeci informacja czy dany snp nie powstal z rozbicia long SNP-a
    vcf_f = vcf.Reader(filename=vcf_file)
    for record in vcf_f:
        if record.CHROM not in lines_from_vcf.keys():
            lines_from_vcf[record.CHROM] = {}
        # We split complex SNPs into "simple SNPs"
        if len(record.REF) == len(record.ALT[0]) and len(record.REF) > 1:
            for i in range(len(record.REF)):
                lines_from_vcf[record.CHROM][record.POS - 1 + i] = [record.REF[i], str(record.ALT[0])[i], "long"]
                # ALT itself is a list even for a simple mutation
        else:
            lines_from_vcf[record.CHROM][record.POS - 1] = [record.REF, str(record.ALT[0]), "short"]
    return lines_from_vcf


if __name__ == '__main__':
    bam_file = sys.argv[1]  # bam file ktorego uzywamy do okreslenia pokrycia i uzycia allelu
    vcf_file = sys.argv[2]  # vcf file, input z mutacjami

    slownik_mutacji = read_vcf(vcf_file=vcf_file)
    slownik_wynikow = {}  # "symulujmey" output freebayes dodatkowo laczymi "long SNPs' w jeden record
    # klucze beda te same nazwa segmentu, a potem pozycja (w notacji VCF-a!), potem lista dwuelementowa z pokryciem
    # i uzyciem alternatywnego wariantu

    for segment in slownik_mutacji.keys():
        in_long_snp = False
        samfile = pysam.AlignmentFile(bam_file, "rb")
        for pileupcolumn in samfile.pileup(segment):
            if pileupcolumn.pos in slownik_mutacji[segment].keys():
                # We create and entry for a give segment in output dict
                if segment not in slownik_wynikow.keys():
                    slownik_wynikow[segment] = {}

                # Here we check if analyzed SNP is a part of a "long SNP"
                # if so we memorize the starting position of this long SNP
                if slownik_mutacji[segment][pileupcolumn.pos][2] == "long" and not in_long_snp:
                    long_snp_position = pileupcolumn.pos
                    in_long_snp = True
                    slownik_wynikow[segment][pileupcolumn.pos + 1] = [[], []]  #create entry for this snp

                elif slownik_mutacji[segment][pileupcolumn.pos][2] == "long" and in_long_snp:
                    in_long_snp = True
                    # we continue long snp we DON'T create a separate entry for this NPS
                elif slownik_mutacji[segment][pileupcolumn.pos][2] == "short":
                    in_long_snp = False
                    slownik_wynikow[segment][pileupcolumn.pos + 1] = [[], []]  # short SNP or indel

                status_mutacji = deterine_status(segment_name=segment,
                                                 position_start=pileupcolumn.pos,
                                                 ref_allel=slownik_mutacji[segment][pileupcolumn.pos][0],
                                                 alt_allel=slownik_mutacji[segment][pileupcolumn.pos][1])

                liczba_ref = 0  # ile odczytow ma ref na tej pozycji dla snp-ow
                liczba_alt = 0  # ile odczytow ma alt na tej pozycji dla snp-ow
                liczba_del = 0  # ile odczytow ma delecje na  KOLEJNEJ pozycji dla indeli
                liczba_ins = 0  # ile odczytow ma insercje na KOLEJNEJ pozycji  dla indeli
                valid = 0  # ilosc odczytow
                for pileupread in pileupcolumn.pileups:
                    valid += 1
                    if not pileupread.is_del:
                        # odczyt nie ma delecji na tej pozycji
                        if status_mutacji == "snp":
                            if (pileupread.alignment.query_sequence[pileupread.query_position] ==
                                    slownik_mutacji[segment][pileupcolumn.pos][0]):
                                liczba_ref += 1
                            elif (pileupread.alignment.query_sequence[pileupread.query_position] ==
                                  slownik_mutacji[segment][pileupcolumn.pos][1]):
                                liczba_alt += 1
                            else:
                                # neither reference nor alt
                                pass

                        if status_mutacji == "ins" or status_mutacji == "del":
                            # For insertion we peek at the next position if the next operation is an insertion, indel will be positive.
                            if pileupread.indel > 0:
                                liczba_ins += 1
                            elif pileupread.indel < 0:
                                liczba_del += 1
                    else:

                        liczba_del += 1

                if in_long_snp:
                        slownik_wynikow[segment][long_snp_position + 1][0].append(valid)
                        slownik_wynikow[segment][long_snp_position + 1][1].append(liczba_alt)
                else:

                    slownik_wynikow[segment][pileupcolumn.pos + 1][0].append(valid)
                    if status_mutacji == "snp":
                        slownik_wynikow[segment][pileupcolumn.pos + 1][1].append(liczba_alt)
                    elif status_mutacji == "ins":
                        slownik_wynikow[segment][pileupcolumn.pos + 1][1].append(liczba_ins)
                    elif status_mutacji == "del":
                        slownik_wynikow[segment][pileupcolumn.pos + 1][1].append(liczba_del)

        samfile.close()
    #  Writing down output form
    #  SEGMENT_NAME\tPOSITION\tDEPTH\tALT_USAGE
    #  print(slownik_wynikow)
    with open('depth_and_usage.txt', "w") as f:
        for segment in slownik_wynikow.keys():
            for position in slownik_wynikow[segment]:
                usage = float(np.mean(slownik_wynikow[segment][position][1]) / np.mean(slownik_wynikow[segment][position][0]))
                f.write(f'{segment}\t{position}\t{int(np.mean(slownik_wynikow[segment][position][0]))}\t{usage:.2f}\n')
