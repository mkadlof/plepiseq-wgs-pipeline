process fastqc {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path('initialfastq_forward_fastqc.txt'), path('initialfastq_reverse_fastqc.txt')

    script:
    """
    fastqc --format fastq \
           --threads ${params.threads} \
           --memory ${params.memory} \
           --extract \
           --delete \
           --outdir . \
           ${reads[0]} ${reads[1]}
    parse_fastqc_output.py ${sampleId}_1_fastqc/fastqc_data.txt ${params.quality_initial} >> initialfastq_forward_fastqc.txt
    parse_fastqc_output.py ${sampleId}_2_fastqc/fastqc_data.txt ${params.quality_initial} >> initialfastq_reverse_fastqc.txt
    """
}