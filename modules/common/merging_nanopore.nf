process merging_nanopore {
    // Ten proces laczy wyniki maskowania dla sciezki z nanopore dla SARS i RSV
    // Gdzie mamy oddzielne mmapowania dla odczytow strict i overshoot
    tag "merging:${sampleId}"
    container  = params.main_image
    input:
    tuple val(sampleId), path('trimmed_first.bam'), path('trimmed_first.bam.bai'),  path(genome), val(QC_status), path('trimmed_second.bam'), path('trimmed_second.bam.bai'),  val(QC_status_2)

    output:
    tuple val(sampleId), path('merged.bam'), path('merged.bam.bai'), path(genome), val(QC_status), emit: to_medaka

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch merged.bam
      touch merged.bam.bai
    else
      samtools merge -o merged_initial.bam trimmed_first.bam trimmed_second.bam
      samtools sort -@ ${params.threads} -o merged.bam merged_initial.bam
      samtools index merged.bam
    fi
    """
}
