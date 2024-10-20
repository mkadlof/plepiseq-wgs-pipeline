process bwa {
    tag "bwa:${sampleId}"
    maxForks 5

    input:
    tuple val(sampleId), path(reads), path(ref_genome_with_index), val(QC_status)

    output:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), env(QC_exit)

    script:
    // Check the index of a file with fasta extension in ref_genome_with_index list
    
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }

    """
    if [ ${QC_status} == "nie" ]; then
      touch mapped_reads.bam
      touch mapped_reads.bam.bai
      QC_exit="nie"
    else
      if [ final_index -eq -1 ]; then
        touch mapped_reads.bam
        touch mapped_reads.bam.bai
        QC_exit="nie" // upstream modules did not provide a valid genome file
      else
        bwa mem -t ${params.threads} -T 30 ${ref_genome_with_index[final_index]} ${reads[0]} ${reads[1]} | \
        samtools view -@ ${params.threads} -Sb -f 3 -F 2048 - | \
        samtools sort -@ ${params.threads} -o mapped_reads.bam -
        samtools index mapped_reads.bam
        NO_READS=`samtools view mapped_reads.bam | wc -l`
        if [ \${NO_READS} -lt 1 ]; then
          QC_exit="nie"
        else
          QC_exit="tak"
        fi  # if na brak poprawnych odczytow po mapowaniu
      fi # if na nie podanie poprawnego genomu
    fi # if na QC status
    """
}
