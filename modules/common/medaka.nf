process medaka_nanopore {
    // Ten proces laczy wyniki maskowania dla sciezki z nanopore dla SARS i RSV
    // Gdzie mamy oddzielne mmapowania dla odczytow strict i overshoot
    tag "merging:${sampleId}"
    container  = params.medaka_image
    input:
    tuple val(sampleId), path('trimmed.bam'), path('trimmed.bam.bai'),  path('genome.fasta'), val(QC_status)

    output:
    tuple val(sampleId), path('merged.bam'), path('merged.bam.bai'), path(genome), val(QC_status), emit: to_medaka

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      # TBD
    else
       medaka consensus --model ${medaka_model} --threads ${cpu} --chunk_len ${chunk_size} --chunk_ovlp ${chunk_overlap} ${prefix}_forvariants.bam ${prefix}_forvariants.hdf >> ${log} 2>&1                                                                                                                                                                                    echo "medaka variant ${input_genome} ${prefix}_forvariants.hdf ${prefix}_medaka.vcf" >> ${log}                          medaka variant ${input_genome} ${prefix}_forvariants.hdf ${prefix}_medaka.vcf >> ${log} 2>&1                                                                                                                                                    #echo "medaka tools annotate ${prefix}_medaka.vcf ${input_genome} ${prefix}_forvariants.bam ${prefix}_medaka_annotated.vcf" >> ${log}                                                                                                           #medaka tools annotate ${prefix}_medaka.vcf ${input_genome} ${prefix}_forvariants.bam ${prefix}_medaka_annotated.vcf >> ${log} 2>&1                                                                                                                                                                                                                                     # podmiana na plik przed maskowaniem primerow do liczenia pokrycia mutacji                                              # bo w samplach 07-09 mamy sytache gdzie delecja jedt na pierwszej pozycji za zamaskowanym primer                       # ,utacja jest wykrywana ale liczac pokrycie delecji program annotate bierze pozycje -1 przed delecja                   # ktora w pliku forvariant jest zamaskowana i ma pokrycie 0                                                             echo "medaka tools annotate ${prefix}_medaka.vcf ${input_genome} ${prefix}_to_clip_sorted.bam ${prefix}_medaka_annotated.vcf" >> ${log}                                                                                                         medaka tools annotate ${prefix}_medaka.vcf ${input_genome} ${prefix}_to_clip_sorted.bam ${prefix}_medaka_annotated.vcf >> ${log} 2>&1

    fi
    """
}
