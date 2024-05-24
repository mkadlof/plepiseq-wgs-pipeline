process coinfection_varscan {
    tag "Varscan during coinfection analysis for sample:\t$sampleId"
    //publishDir "${params.results_dir}//${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(for_contamination_mpileup)

    output:
    tuple val(sampleId), path('detected_variants_varscan_contamination.txt')

    script:
    """
    varscan_qual=`echo "${params.quality_snp} - 1" | bc -l`
    java -jar /opt/varscan/VarScan.v2.4.6.jar pileup2snp ${for_contamination_mpileup} \
                                                     --min-avg-qual \${varscan_qual} \
                                                     --p-value 0.9 \
                                                     --min-var-freq 0.05 \
                                                     --min-coverage ${params.min_cov} \
                                                     --variants \
                                                     --min-reads2 0 > detected_variants_varscan_contamination.txt

    """
}
