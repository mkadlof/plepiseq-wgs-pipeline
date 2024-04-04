process coinfection_ivar {
    publishDir "results/${sampleId}", mode: 'symlink', pattern: 'mapped_reads.bam'

    input:
    tuple val(sampleId), path(mapped_reads),  path(mapped_reads_bai)
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(sampleId), path('forcontaminations_sorted.bam'), path('forcontaminations_sorted.bam.bai'), path('forcontaminations.mpileup')

    script:
    """
    ivar trim -i ${mapped_reads} \
              -b ${params.primers} \
              -m ${params.length} \
              -q ${params.quality_initial} \
              -e \
              -p forcontaminations

    samtools sort -@ ${params.threads} -o forcontaminations_sorted.bam forcontaminations.bam
    samtools index forcontaminations_sorted.bam

    samtools mpileup --max-depth 10000 \
                 --fasta-ref ${reference_fasta} \
                 --min-BQ ${params.quality_SNP} \
                 forcontaminations_sorted.bam >> forcontaminations.mpileup
    """
}