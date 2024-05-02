# Functional analysis: snpEff

```bash
    java -jar /opt/snpEff/snpEff.jar ann -noStats ${params.ref_genome_id} \
         ${consensus_vcf_gz} > detected_variants_consensus_annotated.vcf
```

The analysis was conducted using the snpEFF program. The only required input, besides specifying the genome name present
in the database, is a VCF file containing mutations. This file is then parsed into a text file using bcftools. At this
stage, a consensus VCF file is created. Note that this file is only used for functional analysis, not for creating the
genome sequence.

```Bash
    bcftools query --format '%REF%POS%ALT| %ANN \n' \
               detected_variants_consensus_annotated.vcf.gz | \
                  cut -d "|" -f1,3,5,12 | \
                  tr "|" "\t" | \
                  awk  'BEGIN {OFS = "\t"} {if ( \$2 == "upstream_gene_variant" || \$2 == "downstream_gene_variant") {gene="."; aa="."} else {gene=\$3; aa=\$4}; print gene, \$1, aa, \$2}' > detected_variants_consensus_annotated.txt
```

The `-n+2` option ensures that only mutations occurring in at least two partial files are reported, while `-c` all defines
how positions in different files are treated, with "all" indicating that entries in different files covering the same
position are treated as identical.

Subsequently, the snpEFF program utilizes the files `0000.vcf` and `0001.vcf` created in the dir directory, which contain
VCF files from the lofreq and freebayes programs filtered to include only positions present in at least 2 programs.
These files are parsed in the script to obtain a nicely formatted table.

Note: In the case of complicated variants (deletions of different lengths), there is no guarantee that such a position
will be found in the annotated VCF file, which is a limitation of the isec function.