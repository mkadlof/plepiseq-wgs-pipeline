process viterbi {

    input:
    tuple val(sampleId), path(bam), path(bai)
    path ref_genome

    output:
    tuple val(sampleId), path('forvariants.bam'), path('forvariants.bam.bai')

    script:
    """
    lofreq viterbi --ref ${ref_genome} \
                   --out clean_sort_dedup_trimmed_sort_viterbi.bam \
                   ${bam}

    lofreq indelqual --ref ${ref_genome} \
                     --out forvariants.bam \
                     --dindel clean_sort_dedup_trimmed_sort_viterbi.bam
    samtools sort -@ ${params.threads} \
                  -o forvariants.bam \
                  forvariants.bam
    samtools index forvariants.bam
    """
}