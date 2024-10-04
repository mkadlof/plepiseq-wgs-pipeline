process dehumanization  {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy'

    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), path(reads)

    output:
    tuple val(sampleId), path('forward_paired_nohuman.fq.gz'), path('reverse_paired_nohuman.fq.gz')

    script:
    """
    samtools view mapped_reads.bam | cut -f1 | sort | uniq >> lista_id_nohuman.txt
    seqtk subseq ${reads[0]} lista_id_nohuman.txt >> forward_paired_nohuman.fq
    seqtk subseq ${reads[1]} lista_id_nohuman.txt >> reverse_paired_nohuman.fq

    gzip forward_paired_nohuman.fq
    gzip reverse_paired_nohuman.fq
    """
}
