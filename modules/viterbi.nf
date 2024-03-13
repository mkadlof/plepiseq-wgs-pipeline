process viterbi {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(sampleId), path('forvariants.bam'), path('forvariants.bam.bai')

    script:
    """
    lofreq viterbi --ref ${reference_fasta} \
                   --out clean_sort_dedup_trimmed_sort_viterbi.bam \
                   ${bam}

    lofreq indelqual --ref ${reference_fasta} \
                     --out forvariants.bam \
                     --dindel clean_sort_dedup_trimmed_sort_viterbi.bam
    samtools sort -@ ${params.threads} \
                  -o forvariants.bam \
                  forvariants.bam
    samtools index forvariants.bam
    """
}