process bwa {
    publishDir "results/${sampleId}", mode: 'symlink', pattern: 'mapped_reads.bam'
    publishDir "results/${sampleId}", mode: 'symlink', pattern: 'mapped_reads.bam.bai'

    input:
    tuple val(sampleId), path(reads)
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai')

    script:
    """
    bwa index ${reference_fasta}
    bwa mem -t ${params.threads} -T 30 ${reference_fasta} ${reads[0]} ${reads[1]} | \
        samtools view -@ ${params.threads} -Sb -f 3 -F 2048 - | \
        samtools sort -@ ${params.threads} -o mapped_reads.bam -
    samtools index mapped_reads.bam
    """
}