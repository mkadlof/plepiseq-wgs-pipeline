process freeBayes {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(sampleId), path('detected_variants_freebayes_final.vcf.gz'), path('detected_variants_freebayes_final.vcf.gz.tbi')
    tuple val(sampleId), path('freebayes.fa')

    script:
    """
    freebayes --limit-coverage ${params.max_depth} \
              --min-coverage ${params.min_cov} \
              --min-mapping-quality 20 \
              --min-base-quality ${params.quality_SNP} \
              --use-mapping-quality \
              --fasta-reference ${reference_fasta} \
              --ploidy 1 \
              ${bam} > detected_variants_freebayes.vcf

    cat detected_variants_freebayes.vcf | \
        bcftools norm --check-ref w \
                      --rm-dup all \
                      --fasta-ref ${reference_fasta} | \
                          bcftools norm --check-ref w \
                                        --multiallelics -indels \
                                        --fasta-ref ${reference_fasta} > detected_variants_freebayes_fix.vcf
    
    qual=`echo ${params.pval} | awk '{print int(10*-log(\$1)/log(10))}'`
    
    bcftools filter --include "QUAL >= \${qual} & INFO/DP >= ${params.min_cov} & (SAF  + SAR)/(SRF + SRR + SAF + SAR) > ${params.upper_ambig} " \
             detected_variants_freebayes_fix.vcf > detected_variants_freebayes_fix_high.vcf
    
    bgzip --force detected_variants_freebayes_fix_high.vcf
    tabix detected_variants_freebayes_fix_high.vcf.gz
    
    bcftools filter --include "QUAL >= \${qual} & INFO/DP >=  ${params.min_cov}  & (SAF  + SAR)/(SRF + SRR + SAF + SAR) >= ${params.lower_ambig}  & (SAF  + SAR)/(SRF + SRR + SAF + SAR) <= ${params.upper_ambig} " \
             detected_variants_freebayes_fix.vcf > tmp_low.vcf
    introduce_amb_2_vcf.py tmp_low.vcf \
           detected_variants_freebayes_fix_ambig.vcf

    bgzip --force detected_variants_freebayes_fix_ambig.vcf
    tabix detected_variants_freebayes_fix_ambig.vcf.gz

    bcftools concat detected_variants_freebayes_fix_high.vcf.gz \
                    detected_variants_freebayes_fix_ambig.vcf.gz | \
                        bcftools sort --output-type z > detected_variants_freebayes_final.vcf.gz
    tabix detected_variants_freebayes_final.vcf.gz

    cat ${reference_fasta} | \
            bcftools consensus --samples - \
                               detected_variants_freebayes_final.vcf.gz > freebayes.fa
    """
}