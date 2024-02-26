process consensusAnalysis {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(variants_varscan), path(variants_varscan_tbi),
                         path(variants_freebayes), path(variants_freebayes_tbi),
                         path(variants_lofreq), path(variants_freebayes_tbi)


    output:
    tuple val(sampleId), path('detected_variants_consensus_annotated.txt')

    script:
    """
    bcftools isec --nfiles +2 \
                  --collapse all \
                  --prefix tmp \
                  ${variants_lofreq} ${variants_freebayes} ${variants_varscan}

    mv tmp/* .

    bgzip --force 0000.vcf; tabix 0000.vcf.gz
    bgzip --force 0001.vcf; tabix 0001.vcf.gz

    bcftools merge 0000.vcf.gz \
                   0001.vcf.gz | \
                        bcftools sort --output-type z > detected_variants_consensus_final.vcf.gz

    tabix detected_variants_consensus_final.vcf.gz
    java -jar /opt/snpEff/snpEff.jar ann -noStats ${params.ref_genome_id} \
             detected_variants_consensus_final.vcf.gz > detected_variants_consensus_annotated.vcf

    bgzip --force detected_variants_consensus_annotated.vcf;
    tabix detected_variants_consensus_annotated.vcf.gz

    bcftools query --format '%REF%POS%ALT| %DP | %AF | %ANN \n' \
                   detected_variants_consensus_annotated.vcf.gz | \
                      cut -d "|" -f1,2,3,5,7,14 | \
                      tr "|" "\t" | \
                      awk  'BEGIN {OFS = "\t"} {if ( \$4 == "upstream_gene_variant" || \$4 == "downstream_gene_variant") {gene="."; aa="."} else {gene=\$5; aa=\$6}; print gene, \$1, aa, \$2, \$3*100"%" }' > detected_variants_consensus_annotated.txt

# TODO: The code below has unhandled exception
# In 2.1 pipeline it worked because bash ignored error and moved on, but
# Nextflow fails.

#    i=0
#    lofreq_string=`cat sites.txt | cut -f5 | cut -b1`
#    lofreq_string=(` echo \${lofreq_string} `)
#    freebayes_string=`cat sites.txt | cut -f5 | cut -b2`
#    freebayes_string=(` echo \${freebayes_string} `)
#    varscan_string=`cat sites.txt | cut -f5 | cut -b3`
#    varscan_string=(` echo \${varscan_string} `)
#
#    while read L; do
#        echo -e "\${L}\t\${lofreq_string[\${i}]}\t\${freebayes_string[\${i}]}\t\${varscan_string[\${i}]}" > detected_variants_consensus_annotated.txt.tmp
#        i=\$((i+1))
#    done < detected_variants_consensus_annotated.txt
#
#    cat detected_variants_consensus_annotated.txt.tmp | cut -f1,2,3,6,7,8 > tmp
#    mv tmp detected_variants_consensus_annotated.txt
    """
}