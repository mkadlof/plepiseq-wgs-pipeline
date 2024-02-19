process bwa {
    publishDir "results/${sampleId}", mode: 'symlink', pattern: 'mapped_reads.bam'
    publishDir "results/${sampleId}", mode: 'symlink', pattern: 'mapped_reads.bam.bai'

    input:
    tuple val(sampleId), path(reads)
    path ref_genome

    output:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai')

    script:
    """
    bwa index ${ref_genome}
    bwa mem -t ${params.threads} -T 30 ${ref_genome} ${reads[0]} ${reads[1]} | \
        samtools view -@ ${params.threads} -Sb -f 3 -F 2048 - | \
        samtools sort -@ ${params.threads} -o mapped_reads.bam -
    samtools index mapped_reads.bam
    """
}