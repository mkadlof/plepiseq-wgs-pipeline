process coinfection_ivar {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(mapped_reads),  path(mapped_reads_bai)
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(sampleId), path('for_contamination_sorted.bam'), path('for_contamination_sorted.bam.bai')
    tuple val(sampleId), path('for_contamination.mpileup')

    script:
    """
    ivar trim -i ${mapped_reads} \
              -b ${params.primers} \
              -m ${params.length} \
              -q ${params.quality_initial} \
              -e \
              -p for_contamination

    samtools sort -@ ${params.threads} -o for_contamination_sorted.bam for_contamination.bam
    samtools index for_contamination_sorted.bam

    samtools mpileup --max-depth 10000 \
                 --fasta-ref ${reference_fasta} \
                 --min-BQ ${params.quality_SNP} \
                 for_contamination_sorted.bam >> for_contamination.mpileup
    """
}