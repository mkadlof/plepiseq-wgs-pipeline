process sort_and_index {
    tag "sort_and_index:${sampleId}"

    input:
    tuple val(sampleId), path(bam)

    output:
    tuple val(sampleId), path("${bam}"), path("${bam}.bai")

    script:
    """
    samtools sort -o ${bam} ${bam}
    samtools index ${bam}
    """
}
