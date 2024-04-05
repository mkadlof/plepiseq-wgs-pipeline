#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import sys
import seaborn as sns
from scipy.stats import kstest

path1 = sys.argv[1]  # plik wynikowy z Varscana naszego sample
title = sys.argv[2]
path2 = sys.argv[3]  # sample 9 / zanieczyszczona/
path3 = sys.argv[4]  # sample 17 / zanieczyszczona /
path4 = sys.argv[5]  # sample 32

# przerabiamy wektor procentow z kazdeg
# uzywamy tylko mutacji a nie indeli ...
dane_sample = pd.read_csv(path1, sep='\t')
dane_sample.VarFreq = [float(x.split('%')[0]) / 100 for x in dane_sample.VarFreq]
dane_sample_fortest = [x for x in dane_sample.VarFreq if 0.1 <= x <= 0.9]

# Sample04
dane_sample_04 = pd.read_csv(path2, sep='\t')
dane_sample_04.VarFreq = [float(x.split('%')[0]) / 100 for x in dane_sample_04.VarFreq]
dane_sample_04_fortest = [x for x in dane_sample_04.VarFreq if 0.1 <= x <= 0.9]

# SAmple05
dane_sample_05 = pd.read_csv(path3, sep='\t')
dane_sample_05.VarFreq = [float(x.split('%')[0]) / 100 for x in dane_sample_05.VarFreq]
dane_sample_05_fortest = [x for x in dane_sample_05.VarFreq if 0.1 <= x <= 0.9]

dane_sample_32 = pd.read_csv(path4, sep='\t')
dane_sample_32.VarFreq = [float(x.split('%')[0]) / 100 for x in dane_sample_32.VarFreq]
dane_sample_32_fortest = [x for x in dane_sample_32.VarFreq if 0.1 <= x <= 0.9]

if len(dane_sample_fortest) == 0:
    dane_sample_fortest = [0]

# Rysowanie histogramu z uzyciem allelu alternatywnego
f, ax = plt.subplots(figsize=(7, 5))
sns.despine(f)
sns.histplot(dane_sample_fortest)
ax.set_xlim([0, 1])
ax.set_title(f'Alternative alleles frequencies for {title}')
f.savefig(f'{title}_alternative_alleles_frequencies.png', dpi=600)

# Test ks
pval_sample04 = kstest(dane_sample_fortest, dane_sample_04_fortest).pvalue
pval_sample05 = kstest(dane_sample_fortest, dane_sample_05_fortest).pvalue
pval_sample32 = kstest(dane_sample_fortest, dane_sample_32_fortest).pvalue

# print(pval_sample04, pval_sample05)
# Prosta regula oba p-valu ponizek 0.1 nie zanieczyszczona, oba powyzej zanieczyszczona

with open(f'{title}_contamination_summary.txt', 'w') as f:
    if len(dane_sample_fortest) <= 10:
        f.write(f'{title} is not contaminated with other SARS-CoV-2 material low number\n')
    else:
        if pval_sample04 < 0.1 and pval_sample05 < 0.1 and pval_sample32 < 0.1:
            f.write(f'{title} is not contaminated with other SARS-CoV-2 material\t{pval_sample04}\t{pval_sample05}\t{pval_sample32}\n')
        elif pval_sample04 > 0.1 and pval_sample05 > 0.1 and pval_sample32 > 0.1:
            f.write(f'{title} is contaminated with other SARS-CoV-2 material\t{pval_sample04}\t{pval_sample05}\t{pval_sample32}\n')
        elif pval_sample04 > 0.1 or pval_sample05 > 0.1 or pval_sample32 > 0.1:
            f.write(f'{title} might be contaminated with other SARS-CoV-2 material\t{pval_sample04}\t{pval_sample05}\t{pval_sample32}\n')
