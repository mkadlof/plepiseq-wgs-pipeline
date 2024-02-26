process functionalAnalysis {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(variants_varscan), path(variants_varscan_tbi),
                         path(variants_freebayes), path(variants_freebayes_tbi),
                         path(variants_lofreq), path(variants_freebayes_tbi)

    output:
    tuple val(sampleId), path('detected_variants_varscan_annotated.txt'),
                         path('detected_variants_freebayes_annotated.txt'),
                         path('detected_variants_lofreq_annotated.txt')

    script:
    """
    java -jar /opt/snpEff/snpEff.jar ann -noStats -nodownload -verbose ${params.ref_genome_id} ${variants_varscan} > detected_variants_varscan_annotated.vcf
    java -jar /opt/snpEff/snpEff.jar ann -noStats -nodownload -verbose ${params.ref_genome_id} ${variants_freebayes} > detected_variants_freebayes_annotated.vcf
    java -jar /opt/snpEff/snpEff.jar ann -noStats -nodownload -verbose ${params.ref_genome_id} ${variants_lofreq} > detected_variants_lofreq_annotated.vcf
    
    bgzip --force detected_variants_varscan_annotated.vcf
    tabix detected_variants_varscan_annotated.vcf.gz
    bgzip --force detected_variants_freebayes_annotated.vcf
    tabix detected_variants_freebayes_annotated.vcf.gz
    bgzip --force detected_variants_lofreq_annotated.vcf
    tabix detected_variants_lofreq_annotated.vcf.gz
    
    bcftools query -f '%REF%POS%ALT| %DP | %AF | %ANN \n' detected_variants_varscan_annotated.vcf.gz | \
        cut -d "|" -f1,2,3,5,7,14 | \
        tr "|" "\t" | \
        awk 'BEGIN {OFS = "\t"} {if ( \$4 == "upstream_gene_variant" || \$4 == "downstream_gene_variant") {gene="."; aa="."} else {gene=\$5; aa=\$6}; print gene, \$1, aa, \$2, \$3*100"%" }' > detected_variants_varscan_annotated.txt
    bcftools query -f '%REF%POS%ALT | %DP |  %SAF |  %SAR | %SRF | %SRR |  %ANN \n' detected_variants_freebayes_annotated.vcf.gz | \
        cut -d "|" -f1,2,3,4,5,6,8,10,17 | \
        tr "|" "\t" | \
        awk  'BEGIN {OFS = "\t"}  {if ( \$7 == "upstream_gene_variant" || \$7 == "downstream_gene_variant") {gene="."; aa="."} else {gene=\$8; aa=\$9}; print gene, \$1, aa, \$2, (\$3+\$4)/(\$3+\$4+\$5+\$6)*100"%" }' > detected_variants_freebayes_annotated.txt
    ls detected_variants_lofreq_annotated.vcf.gz
    bcftools query -f '%REF%POS%ALT| %DP | %AF | %ANN \n' detected_variants_lofreq_annotated.vcf.gz | \
        cut -d "|" -f1,2,3,5,7,14 | \
        tr "|" "\t" | \
        awk 'BEGIN {OFS = "\t"} {if ( \$4 == "upstream_gene_variant" || \$4 == "downstream_gene_variant") {gene="."; aa="."} else {gene=\$5; aa=\$6}; print gene, \$1, aa, \$2, \$3*100"%" }' > detected_variants_lofreq_annotated.txt
    """
}