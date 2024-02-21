process wgsMetrics {

    input:
    tuple val(sampleId), path(bam), path(bai)
    path ref_genome

    output:
    tuple val(sampleId), path('picard_statistics.txt')

    script:
    """
    java -jar /opt/picard/picard.jar CollectWgsMetrics --REFERENCE_SEQUENCE ${ref_genome} \
                                                   --MINIMUM_BASE_QUALITY ${params.quality_initial} \
                                                   --MINIMUM_MAPPING_QUALITY 30 \
                                                   --INPUT ${bam} \
                                                   --OUTPUT picard_statistics.txt
    """
}