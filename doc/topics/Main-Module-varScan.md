# Main: Module varScan

Below is the procedure on how we use and parse the varscan program:
```Bash
    samtools mpileup --max-depth ${params.max_depth} \
                 --fasta-ref ${reference_fasta} \
                 --min-BQ ${params.quality_snp} \
                 ${bam} >> ${bam}.mpileup
```
Creating the mpileup file required by the program. The option `-d` (note: properly functioning from `samtools` version 1.09 onwards) ignores further reads mapping to a given position if the coverage exceeds the assumed value.

```Bash
    varscan_qual=`echo "${params.quality_snp} - 1" | bc -l`
    java -jar /opt/varscan/VarScan.v2.4.6.jar pileup2cns ${bam}.mpileup \
         --min-avg-qual \${varscan_qual} \
         --p-value ${params.pval} \
         --min-var-freq ${params.lower_ambig} \
         --min-coverage ${params.min_cov} \
         --variants \
         --min-reads2 0 > detected_variants_varscan.txt
```

Calling varscan, the above options signify:

 - `--min-coverage` – minimum coverage at a given position required for the program to identify a variant.
 - `--min-reads2` - minimum number of reads supporting the variant, set to 0.
 - `--min-avg-qual` – minimum read quality at this position required to include the read in the analysis.
 - `--p-value` – threshold p-value the variant must achieve to be reported.
 - `--min-var-freq` – minimum percentage of reads with a non-reference allele required to report a variant.
 - `--variants` – besides SNPs, the program should also identify short INDELs.`

```Bash
    parse_vcf_output_final.py detected_variants_varscan.txt ${params.upper_ambig} ${params.pval}
```

A parser converting the text file `detected_variants_varscan.txt` to a vcf file format. Compared to the txt file, the vcf file introduces the following changes:

1. Positions reported as deletions/insertions relative to the reference genome, e.g.,
    
       MN908947.3    22204    T    +GAGCCAGAA

   Containing symbols "+" or "-" are corrected to:

       MN908947.3    22204    .    T    TGAGCCAGAA

2. Heterozygous positions, e.g.,

       MN908947.3    21766    .    ACATGTC    A

3. Positions where the frequency of using the alternative allele is greater than the value provided in the variable $upper_ambig are converted from ambiguous to the alternative allele version. Varscan returns ambiguous positions even with a reference allele frequency of 0.75, e.g.,
    
       MN908947.    25324    C    M (where the frequency of allele A is 60%)
    
   are corrected to

       MN908947.3    25324    .    C    A

4. The QUAL field, which does not have a clear representation in the varscan output, is changed to 30 in the vcf file.

```Bash
bcftools norm --check-ref w \
                  --rm-dup all \
                  --fasta-ref ${reference_fasta}\
                   detected_variants_varscan.vcf.gz | \
                       bcftools norm --check-ref w \
                                     --multiallelics -indels \
                                     --fasta-ref ${reference_fasta} | \
                                           bcftools filter \
                                                    --include "QUAL >= \${qual} && AF >= ${params.lower_ambig} && DP >= ${params.min_cov}" > detected_variants_varscan_final.vcf
```

"Sorting" the vcf file using the norm function from the bcftools package. The `-c w` option causes the program to return a warning instead of an error if ambiguous positions are present in the file, and further process the vcf file. The "-d all" removes duplicates, the `-m -indels` option splits multiallelic positions into individual entries. Generally, vcf files show the richness of changes identified by programs. Sometimes overlapping changes occur, sometimes programs return unusual mutation notations, multiallelic positions, etc. To organize such a vcf file, we call the above command. Then we ensure that the resulting records still meet the quality and coverage criteria. In the case of varscan, this is not as important because when parsing its output, we set the quality to 1 or 30. Simple examples of changes identified in sample 1 from EQA for the freebayes program.

Before:

    MN908947.3    22204    .    TGA    TGAGCCAGAAGA

After:

    MN908947.3    22204    .    T    TGAGCCAGAA

Or splitting a complex mutation
Before:

    MN908947.3    25334    .    GAAG    CACG,CACC,GACC,GAAC

After, splitting the mutation into components

    MN908947.3    25334    .    GAA    CAC
    MN908947.3    25334    .    GAAG    CACC
    MN908947.3    25336    .    AG    CC
    MN908947.3    25337    .    G    C

```Bash
    cat ${reference_fasta} | bcftools consensus --samples - detected_variants_varscan_final.vcf.gz > varscan.fa
```

The above command incorporates all mutations present in the vcf file into the reference genome. The `-s -` option means that we ignore the genotype for the sample (we do not have multiple samples in the vcf file), and we introduce all mutations present in the vcf file. The result is a fasta file format with the genome.