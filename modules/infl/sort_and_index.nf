process sort_and_index {
    tag "sort_and_index:${sampleId}"

    input:
    tuple val(sampleId), path(bam)

    output:
    tuple val(sampleId), path("${bam.baseName}_sorted.bam"), path("${bam.baseName}_sorted.bam.bai")

    script:
    def newBam = "${bam.baseName}_sorted.bam"
    """
    samtools sort -o ${newBam} ${bam}
    samtools index ${newBam}
    """
}
