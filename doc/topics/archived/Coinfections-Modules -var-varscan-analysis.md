# Coinfections: Modules ivar, varscan, analysis

The co-infection section consists of a preparatory module whose task is to trim primer sequences from the sequences.

```Bash
    ivar trim -i ${mapped_reads} \
              -b ${primers} \
              -m ${params.length} \
              -q ${params.quality_initial} \
              -e \
              -p for_contamination
    (...)
    samtools mpileup --max-depth 10000 \
                 --fasta-ref ${reference_fasta} \
                 --min-BQ ${params.quality_snp} \
                 for_contamination_sorted.bam >> for_contamination.mpileup
```

Next, the result is passed to the VarScan program, which identifies mutations in unfiltered samples.

```Bash
    varscan_qual=`echo "${params.quality_snp} - 1" | bc -l`
    java -jar /opt/varscan/VarScan.v2.4.6.jar pileup2snp ${for_contamination_mpileup} \
                                                     --min-avg-qual \${varscan_qual} \
                                                     --p-value 0.9 \
                                                     --min-var-freq 0.05 \
                                                     --min-coverage 20 \
                                                     --variants \
                                                     --min-reads2 0 > detected_variants_varscan_contamination.txt

```

The final module compares allele frequency distributions with known samples that have been identified as co-infected (based on EQA23 tests) and ultimately answers the question about the probability of co-infection occurrence.

