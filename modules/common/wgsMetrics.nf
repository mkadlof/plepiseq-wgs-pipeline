process wgsMetrics {
    tag "wgsMetrics:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "picard_statistics.txt"

    input:
    tuple val(sampleId), path(bam), path(bai), path(ref_genome)

    output:
    tuple val(sampleId), path('picard_statistics.txt')
    tuple val(sampleId), path("wgsMetrics.json")

    script:
    """
    java -jar /opt/picard/picard.jar CollectWgsMetrics --REFERENCE_SEQUENCE ${ref_genome} \
                                                   --MINIMUM_BASE_QUALITY ${params.quality_initial} \
                                                   --MINIMUM_MAPPING_QUALITY ${params.min_mapq} \
                                                   --INPUT ${bam} \
                                                   --OUTPUT picard_statistics.txt

    parse_wgsMetrics.py picard_statistics.txt wgsMetrics.json
    """
}
