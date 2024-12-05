#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import kstest

path1 = sys.argv[1]  # data for an analyzed sample
title = sys.argv[2]  # plot title
list_of_coinfected_samples = sys.argv[3:]  # list of TXT files obtained with Varscan on a BAM file for known coinfected
# samples  prior to any filtering, but with primer sites masked with ivar


# Tested samples data, we extract "VarFreq" column from a file and convert these values to floats
# only allels with 0.1 to 0.9 usage are kept
data_sample = pd.read_csv(path1, sep='\t')
data_sample.VarFreq = [float(str(x).strip('%')) / 100 for x in data_sample.VarFreq]
data_sample_fortest = [x for x in data_sample.VarFreq if 0.1 <= x <= 0.9]

# plot histogram for json output
with open('allele_usage_histogram.txt', 'w') as f:
    f.write(f'#indeks,allel_uzycie,liczbnosc\n')
    i = 1
    ilosc, zakres = np.histogram(data_sample_fortest,
                                 bins=16,
                                 range=(0.1, 0.9),
                                 density=False)
    for a, b in zip(ilosc, zakres):
        f.write(f'{i},{b:.2f},{a}\n')
        i += 1


# if data_sample_fortest is empty (no sites which alternative allele usage is between 0.1 to 0.9) add a dummy value

if len(data_sample_fortest) == 0:
    data_sample_fortest = [0]

# A simple plot for visualization
f, ax = plt.subplots(figsize=(7, 5))
sns.despine(f)
sns.histplot(data_sample_fortest)
ax.set_xlim([0, 1])
ax.set_title(f'Alternative alleles frequencies for {title}')
f.savefig(f'{title}_alternative_alleles_frequencies.png', dpi=600)

# Compare an analyzed sample with known coinfected samples, calculate p-value using KS test
pval_list = []
for known_sample in list_of_coinfected_samples:
    coinfected_sample = pd.read_csv(known_sample, sep='\t')
    coinfected_sample.VarFreq = [float(str(x).strip('%')) / 100 for x in coinfected_sample.VarFreq]
    coinfected_sample_fortest = [x for x in coinfected_sample.VarFreq if 0.1 <= x <= 0.9]

    # run ks test to determine if distributions
    pval = kstest(data_sample_fortest, coinfected_sample_fortest).pvalue
    pval_list.append(pval)

pval_list = np.array(pval_list)
text = "\t".join(map(str, pval_list))

with open(f'{title}_coinfection_summary.txt', 'w') as f:
    if len(data_sample_fortest) <= 10:
        f.write(f'{title} is not coinfected with other SARS-CoV-2 material, less than 10 mixed sites \n')
        print(f"nie 0")
    else:
        if np.all(pval_list < 0.1):
            print(f"nie {np.max(pval_list)}")
            f.write(f'{title} is not coinfected with other SARS-CoV-2 material. '
                    f'P-values to known co-infected samples are {text}\n')
        elif np.all(pval_list > 0.1):
            print(f"tak {np.max(pval_list)}")
            f.write(f'{title} is coinfected with other SARS-CoV-2 material. '
                    f'P-values to known co-infected samples are {text}\n')
        elif np.any(pval_list > 0.1):
            print(f"prawdopodobny {np.max(pval_list)}")
            f.write(f'{title} might be coinfected with other SARS-CoV-2 material it shows similarity to at least one '
                    f'sample that was previously classified as co-infected. P-values to known co-infected '
                    f'samples are {text}\n')
