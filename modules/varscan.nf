process varScan {
    publishDir "results/variants", mode: 'symlink'

    input:
    tuple val(sampleId), path(bam), path(bai), path(lowcoverage_masked_fa)
    path ref_genome

    output:
    tuple val(sampleId), path('detected_variants_varscan.txt'), path('varscan_masked.fa')

    script:
    """
    samtools mpileup --max-depth ${params.max_depth} \
                 --fasta-ref ${ref_genome} \
                 --min-BQ ${params.quality_snp} \
                 ${bam} >> ${bam}.mpileup

    varscan_qual=`echo "${params.quality_snp} - 1" | bc -l`
    java -jar /opt/varscan/VarScan.v2.4.6.jar pileup2cns ${bam}.mpileup \
                                                         --min-avg-qual \${varscan_qual} \
                                                         --p-value ${params.pval} \
                                                         --min-var-freq ${params.lower_ambig} \
                                                         --min-coverage ${params.min_cov} \
                                                         --variants \
                                                         --min-reads2 0 > detected_variants_varscan.txt

    parse_vcf_output_final.py detected_variants_varscan.txt ${params.upper_ambig} ${params.pval}

    bgzip --force detected_variants_varscan.vcf
    tabix detected_variants_varscan.vcf.gz

    qual=`echo ${params.pval} | awk '{print int(10*-log(\$1)/log(10))}'`

    bcftools norm --check-ref w \
                  --rm-dup all \
                  --fasta-ref ${ref_genome}\
                   detected_variants_varscan.vcf.gz | \
                       bcftools norm --check-ref w \
                                     --multiallelics -indels \
                                     --fasta-ref ${ref_genome} | \
                                           bcftools filter \
                                                    --include "QUAL >= \${qual} && AF >= ${params.lower_ambig} && DP >= ${params.min_cov}" > detected_variants_varscan_final.vcf

    bgzip --force detected_variants_varscan_final.vcf
    tabix detected_variants_varscan_final.vcf.gz

    cat ${ref_genome} | bcftools consensus --samples - detected_variants_varscan_final.vcf.gz > varscan.fa
    cat ${lowcoverage_masked_fa} varscan.fa >> tmp_varscan.fa
    mafft --auto --inputorder --quiet tmp_varscan.fa >> tmp_varscan_aln.fa

    get_N.py tmp_varscan_aln.fa
    mv output_varscan_masked.fa varscan_masked.fa
    """
}