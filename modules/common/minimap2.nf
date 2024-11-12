process minimap2 {
    container  = params.main_image 
    tag "minimap2:${sampleId}"
    maxForks 5

    input:
    tuple val(sampleId), path(reads), path("ref_genome.fasta"), path("primers.bed"), val(QC_status)

    output:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), env(QC_exit), emit: only_bam
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), path("ref_genome.fasta"), env(QC_exit), emit: bam_and_genome
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), path("ref_genome.fasta"), path("primers.bed"), env(QC_exit), emit: bam_and_genome_and_primers


    script:

    """
    if [ ${QC_status} == "nie" ]; then
      touch mapped_reads.bam
      touch mapped_reads.bam.bai
      QC_exit="nie"
    else
      minimap2 -a -x map-ont -t ${params.threads} -o tmp.sam ref_genome.fasta ${reads}
      samtools view -@ ${params.threads} -Sb -F 2052 tmp.sam | \
      samtools sort -@ ${params.threads} -o mapped_reads.bam -
      samtools index mapped_reads.bam
      rm tmp.sam
      NO_READS=`samtools view mapped_reads.bam | wc -l`
      if [ \${NO_READS} -lt ${params.min_number_of_reads} ]; then
        QC_exit="nie"
      else
        QC_exit="tak"
      fi  # if na brak poprawnych odczytow po mapowaniu
    fi # if na QC status
    """
}
