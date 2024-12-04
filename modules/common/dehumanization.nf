process dehumanization_illumina  {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "*nohuman.fastq.gz"
    container  = params.main_image

    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), val(QC_status), path(reads)

    output:
    tuple val(sampleId), path("*nohuman.fastq.gz"), emit: to_pubdir
    tuple val(sampleId), path('list_of_dehumanzed_fastas.txt'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch ${sampleId}_forward_paired_nohuman.fastq.gz ${sampleId}_reverse_paired_nohuman.fastq.gz
    else
        samtools view mapped_reads.bam | cut -f1 | sort | uniq > lista_id_nohuman.txt
        seqtk subseq ${reads[0]} lista_id_nohuman.txt | gzip > ${sampleId}_forward_paired_nohuman.fastq.gz
        seqtk subseq ${reads[1]} lista_id_nohuman.txt | gzip > ${sampleId}_reverse_paired_nohuman.fastq.gz
    fi
    for K in `ls *paired_nohuman*gz`; do
        echo -e "${params.results_dir}/${sampleId}/\${K}" >> list_of_dehumanzed_fastas.txt
    done
    """
}

process dehumanization_nanopore {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "*nohuman.fastq.gz"
    container  = params.main_image
    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), val(QC_status), path(reads)

    output:
    tuple val(sampleId), path("*nohuman.fastq.gz"), emit: to_pubdir
    tuple val(sampleId), path('list_of_dehumanzed_fastas.txt'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch ${sampleId}_nohuman.fastq.gz
    else
      samtools view mapped_reads.bam | cut -f1 | sort | uniq >> lista_id_nohuman.txt
      seqtk subseq ${reads} lista_id_nohuman.txt | gzip >> ${sampleId}_nohuman.fastq.gz
    fi
    for K in `ls *nohuman*gz`; do
        echo -e "${params.results_dir}/${sampleId}/\${K}" >> list_of_dehumanzed_fastas.txt
    done
    """
}
