process wgsMetrics {
    tag "wgsMetrics:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "picard_statistics.txt"

    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple val(sampleId2), path(ref_genome)

    output:
    tuple val(sampleId), path('picard_statistics.txt')

    script:
    """
    java -jar /opt/picard/picard.jar CollectWgsMetrics --REFERENCE_SEQUENCE ${ref_genome} \
                                                   --MINIMUM_BASE_QUALITY ${params.quality_initial} \
                                                   --MINIMUM_MAPPING_QUALITY ${params.min_mapq} \
                                                   --INPUT ${bam} \
                                                   --OUTPUT picard_statistics.txt
    """
}
