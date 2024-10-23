process dehumanization  {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy'

    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), path(reads)

    output:
    tuple val(sampleId), path('forward_paired_nohuman.fastq.gz'), path('reverse_paired_nohuman.fastq.gz')

    script:
    """
    samtools view mapped_reads.bam | cut -f1 | sort | uniq > lista_id_nohuman.txt
    seqtk subseq ${reads[0]} lista_id_nohuman.txt | gzip > forward_paired_nohuman.fastq.gz
    seqtk subseq ${reads[1]} lista_id_nohuman.txt | gzip > reverse_paired_nohuman.fastq.gz
    """
}
