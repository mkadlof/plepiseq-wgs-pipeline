process merging {
    tag "Merging alignments for sample:\t$sampleId"
    // publishDir "${params.results_dir}/${sampleId}/filtered_bam", mode: 'symlink', pattern "clean_sort_dedup_trimmed_sort.bam*"

    input:
    tuple val(sampleId), path(filtering_bam), path(ivar_bam)
    path primers
    path pairs

    output:
    tuple val(sampleId), path('clean_sort_dedup_trimmed_sort.bam'), path('clean_sort_dedup_trimmed_sort.bam.bai')

    script:
    """
    samtools merge -o clean_sort_dedup_trimmed_sort_tmp.bam ${filtering_bam} ${ivar_bam}
    samtools sort -@ ${params.threads} -o clean_sort_dedup_trimmed_sort.bam clean_sort_dedup_trimmed_sort_tmp.bam
    samtools index clean_sort_dedup_trimmed_sort.bam
    """
}
