# Main: Module freeBayes

```Bash
    freebayes --limit-coverage ${params.max_depth} \
              --min-coverage ${params.min_cov} \
              --min-mapping-quality 20 \
              --min-base-quality ${params.quality_snp} \
              --use-mapping-quality \
              --fasta-reference ${reference_fasta} \
              --ploidy 1 \
              ${bam} > detected_variants_freebayes.vcf
```

Calling the freebayes program. The options used are:

 - `--limit-coverage ${max_depth}` – limit coverage at a position to this value when identifying variants. Please note that quasi-downsampling used in Step 4 for read filtering does not guarantee that a given position will not exceed the value specified in the max_depth variable.
 - `--min-coverage ${min_cov}` – minimum coverage for a position to be considered as a variant.
 - `-m 20` – alignment quality of a read for its nucleotides to be considered in variant analysis.
 - `-q ${quality_SNP}` – minimum nucleotide quality at a position in a read to be considered in variant counting.
 - `-p 1` – ploidy of the analyzed organism.
 - `-f ${input_genome}` – path to the reference genome.
 - `-j` – use mapping qualities when calculating variant likelihood.

```Bash
    cat detected_variants_freebayes.vcf | \
        bcftools norm --check-ref w \
                      --rm-dup all \
                      --fasta-ref ${reference_fasta} | \
                          bcftools norm --check-ref w \
                                        --multiallelics -indels \
                                        --fasta-ref ${reference_fasta} > detected_variants_freebayes_fix.vcf
```

Normalization of the VCF file similar to [varScan](Main-Module-varScan.md).

```Bash
    bcftools filter --include "QUAL >= \${qual} & INFO/DP >= ${params.min_cov} & (SAF  + SAR)/(SRF + SRR + SAF + SAR) > ${params.upper_ambig} " \
             detected_variants_freebayes_fix.vcf > detected_variants_freebayes_fix_high.vcf
```
   
Freebayes cannot identify ambiguous positions, and at the position corresponding to the alternative allele in the VCF file, it always introduces a symbol: `A`, `T`, `G`, or `C`. Introducing the allele ambiguous symbol at such a position must be done manually. In the first step, we save to the `detected_variants_freebayes_fix_high.vcf` file all mutations in which the frequency of using the alternative allele exceeds the value provided in the `${upper_ambig}` variable and whose quality is greater than 15 (meaning p-value is less than 0.03).

```Bash
    bcftools filter --include "QUAL >= \${qual} & INFO/DP >=  ${params.min_cov}  & (SAF  + SAR)/(SRF + SRR + SAF + SAR) >= ${params.lower_ambig}  & (SAF  + SAR)/(SRF + SRR + SAF + SAR) <= ${params.upper_ambig} " \
             detected_variants_freebayes_fix.vcf > tmp_low.vcf
    introduce_amb_2_vcf.py tmp_low.vcf \
           detected_variants_freebayes_fix_ambig.vcf
```

Next, we extract mutations with allele alternative frequency between the values specified in the variables `$lower_ambig` and `$upper_ambig`. For these positions, we apply a simple Python script where instead of introducing the alt allele symbol, we introduce the ambiguous nucleotide symbol.

```Bash
    bcftools concat detected_variants_freebayes_fix_high.vcf.gz \
                    detected_variants_freebayes_fix_ambig.vcf.gz | \
                        bcftools sort --output-type z > detected_variants_freebayes_final.vcf.gz
```


