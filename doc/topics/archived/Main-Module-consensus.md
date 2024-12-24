# Main: Module consensus

```Bash
    cat ${freebayes_fa} ${lofreq_fa} ${varscan_fa} > tmp_to_consensus.fa
    mafft --auto --inputorder --quiet tmp_to_consensus.fa > tmp_to_consensus_aln.fa
    gen_consensus_seq_final.py tmp_to_consensus_aln.fa
```

Generating consensus from three programs utilizes the gap_consensus function from the Biopython package called by the `gen_consensus_seq_final.py` script. In this function, the threshold value is set to 0.6. This means that at a given position in the alignment, the consensus is considered to be the nucleotide that occurs in at least 60% of the sequences. If there is no nucleotide meeting this criterion, "X" is entered. The Biopython function understands deletions, and they are also recorded in the consensus sequence. Since we have three sequences at each position, we can get values of 0.33 (each sequence at this position has a different allele, which is possible when one program returns the reference allele, another program the alternative allele, and the third program ambiguous positions), 0.66 (2 out of 3 sequences share a common allele), or 1 (all sequences have the same allele). Thus, a threshold of 0.6 means that each position is a consensus of at least two, any two, programs. The resulting consensus sequence then undergoes the standard procedure of masking regions with low coverage.
