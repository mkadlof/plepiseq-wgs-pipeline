process fastqc {
    tag "${prefix} QC analysis for sample:\t$sampleId"
    publishDir "${params.results_dir}/${sampleId}/QC", mode: 'symlink'

    input:
    tuple val(sampleId), path(reads)
    val(prefix)

    output:
    tuple val(sampleId), path("${prefix}_forward_fastqc.txt"), path("${prefix}_reverse_fastqc.txt")

    script:
    """
    fastqc --format fastq \
           --threads ${params.threads} \
           --memory ${params.memory} \
           --extract \
           --delete \
           --outdir . \
           ${reads[0]} ${reads[1]}
    r1=\$(basename ${reads[0]} .fastq.gz)
    r2=\$(basename ${reads[1]} .fastq.gz)
    parse_fastqc_output.py \${r1}_fastqc/fastqc_data.txt ${params.quality_initial} >> ${prefix}_forward_fastqc.txt
    parse_fastqc_output.py \${r2}_fastqc/fastqc_data.txt ${params.quality_initial} >> ${prefix}_reverse_fastqc.txt
    """
}
