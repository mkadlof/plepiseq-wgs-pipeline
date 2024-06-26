process coinfection_ivar {
    tag "coinfection_ivar:${sampleId}"

    input:
    tuple val(sampleId), path(mapped_reads), path(mapped_reads_bai)
    tuple val(sampleId2), path(ref_genome)
    tuple val(sampleId3), path(primers)

    output:
    tuple val(sampleId), path('for_contamination_sorted.bam'), path('for_contamination_sorted.bam.bai')
    tuple val(sampleId), path('for_contamination.mpileup')

    script:
    """
    ivar trim -i ${mapped_reads} \
              -b ${primers} \
              -m ${params.length} \
              -q ${params.quality_initial} \
              -e \
              -p for_contamination

    samtools sort -@ ${params.threads} -o for_contamination_sorted.bam for_contamination.bam
    samtools index for_contamination_sorted.bam
    samtools mpileup -B --max-depth 20000 \
                 --fasta-ref ${ref_genome} \
                 --min-BQ ${params.quality_snp} \
                 for_contamination_sorted.bam >> for_contamination.mpileup
    """
}
