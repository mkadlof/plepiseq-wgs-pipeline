process dehumanization_illumina  {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}*nohuman.fastq.gz"
    container  = params.main_image

    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), val(QC_status), path(reads)

    output:
    tuple val(sampleId), path("*nohuman.fastq.gz"), emit: to_pubdir
    tuple val(sampleId), path('dehumanzed.json'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch nohuman.fastq.gz
      ERR_MSG="This module was eneterd with failed QC and poduced no valid output"
      parse_dehumanization.py --status "nie" --error "\${ERR_MSG}" -o dehumanzed.json
    else
        samtools view mapped_reads.bam | cut -f1 | sort | uniq > lista_id_nohuman.txt
        seqtk subseq ${reads[0]} lista_id_nohuman.txt | gzip > ${sampleId}_forward_paired_nohuman.fastq.gz
        seqtk subseq ${reads[1]} lista_id_nohuman.txt | gzip > ${sampleId}_reverse_paired_nohuman.fastq.gz
        
        find . -name "*paired_nohuman*gz" >> list_of_dehumanzed_fastas.txt
        parse_dehumanization.py --status "tak" \
                                --input_fastas_list list_of_dehumanzed_fastas.txt \
                                --output_path "${params.results_dir}/${sampleId}" \
                                --output dehumanzed.json

    fi
    """
}

process dehumanization_nanopore {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}*nohuman.fastq.gz"
    container  = params.main_image
    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), val(QC_status), path(reads)

    output:
    tuple val(sampleId), path("*nohuman.fastq.gz"), emit: to_pubdir
    tuple val(sampleId), path('dehumanzed.json'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch nohuman.fastq.gz
      ERR_MSG="This module was eneterd with failed QC and poduced no valid output"
      parse_dehumanization.py --status "nie" --error "\${ERR_MSG}" -o dehumanzed.json
    else
      samtools view mapped_reads.bam | cut -f1 | sort | uniq >> lista_id_nohuman.txt
      seqtk subseq ${reads} lista_id_nohuman.txt | gzip >> ${sampleId}_nohuman.fastq.gz
      find . -name "*nohuman.fastq.gz" >> list_of_dehumanzed_fastas.txt
     parse_dehumanization.py --status "tak" \
                             --input_fastas_list list_of_dehumanzed_fastas.txt \
                             --output_path "${params.results_dir}/${sampleId}" \
                             --output dehumanzed.json

    fi
    """
}
