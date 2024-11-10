process medaka {
    // Ten proces laczy wyniki maskowania dla sciezki z nanopore dla SARS i RSV
    // Gdzie mamy oddzielne mmapowania dla odczytow strict i overshoot
    tag "merging:${sampleId}"
    container  = params.medaka_image
    input:
    tuple val(sampleId), path('trimmed.bam'), path('trimmed.bam.bai'),  path('genome.fasta'), val(QC_status)
    val(round)
    output:
    tuple val(sampleId), path('medaka_annotated_filtered.vcf.gz'), path('medaka_annotated_filtered.vcf.gz.tbi'), val(QC_status), emit: vcf

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch medaka_annotated_filtered.vcf.gz
      touch medaka_annotated_filtered.vcf.gz.tbi
    else
      medaka inference --model ${params.medaka_model} \
                       --threads ${params.threads} \
                       --chunk_len ${params.medaka_chunk_len} \
                       --chunk_ovlp ${params.medaka_chunk_overlap} \
                       trimmed.bam forvariants.hdf
     
      medaka vcf forvariants.hdf genome.fasta medaka_initial.vcf
      medaka tools annotate medaka_initial.vcf genome.fasta trimmed.bam medaka_annotated.vcf

      bgzip medaka_annotated.vcf; tabix medaka_annotated.vcf.gz
      if [ ${round} -eq 1 ]; then
        QUAL=`echo ${params.first_round_pval} | awk '{print 10*-log(\$1)/log(10)}'` # dodanie int powodowalo blad w zaookraglaniu, teraz dla 0.1 program poprawnie zwoci wartosc 10
      else
        QUAL=`echo ${params.second_round_pval} | awk '{print 10*-log(\$1)/log(10)}'`
      fi

      bcftools filter -O z -o medaka_annotated_filtered.vcf.gz -i "GQ >= \${QUAL} && DP >= ${params.min_cov}" medaka_annotated.vcf.gz
      tabix medaka_annotated_filtered.vcf.gz
    fi
    """
}
