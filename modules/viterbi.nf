process viterbi {

    input:
    tuple val(sampleId), path(bam), path(bai)
    path ref_genome

    output:
    tuple val(sampleId), path('downsample.bam'), path('downsample.bam.bai')

    script:
    """
    lofreq viterbi --ref ${ref_genome} \
                   --out clean_sort_dedup_trimmed_sort_viterbi.bam \
                   ${bam}

#                    ${OUTPUT_DIR}/${prefix}_clean_sort_dedup_trimmed_sort.bam

    lofreq indelqual --ref ${REFERENCE_GENOME_FASTA} \
                     --out ${LOFREQ_OUTPUT_DIR}/${prefix}_clean_sort_dedup_trimmed_sort_viterbi_lofreq.bam \
                     --dindel ${LOFREQ_OUTPUT_DIR}/${prefix}_clean_sort_dedup_trimmed_sort_viterbi.bam
    samtools sort -@ ${cpu} \
                  -o ${OUTPUT_DIR}/${prefix}_forvariants.bam \
                  ${LOFREQ_OUTPUT_DIR}/${prefix}_clean_sort_dedup_trimmed_sort_viterbi_lofreq.bam
    samtools index ${OUTPUT_DIR}/${prefix}_forvariants.bam
    """
}