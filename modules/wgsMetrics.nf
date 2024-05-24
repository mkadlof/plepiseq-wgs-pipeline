process wgsMetrics {
    tag "Simple statistics for final bam file for sample:\t$sampleId"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "picard_statistics.txt"
    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(sampleId), path('picard_statistics.txt')

    script:
    """
    java -jar /opt/picard/picard.jar CollectWgsMetrics --REFERENCE_SEQUENCE ${reference_fasta} \
                                                   --MINIMUM_BASE_QUALITY ${params.quality_initial} \
                                                   --MINIMUM_MAPPING_QUALITY ${params.mapping_quality} \
                                                   --INPUT ${bam} \
                                                   --OUTPUT picard_statistics.txt
    """
}
