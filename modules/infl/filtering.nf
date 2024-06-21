process filtering {
    tag "filtering:${sampleId}"

    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple val(sampleId), path(ref_genome_fasta)
    tuple val(sampleId), path(primers)

    output:
//     tuple val(sampleId), path('first_pass_sorted.bam'), path('first_pass_sorted.bam.bai')
//     tuple val(sampleId), path('two_amplicons_sorted.bam')
//     tuple val(sampleId), path('Statystyki.txt')

    script:
    """
    # Custom filtering for influenza
    # The script takes in sequence:
    # - BAM for filtering and downsampling
    # - Primer scheme
    # - Target coverage per segment (integer)
    # - Minimum read mapping quality (integer)
    # - Reference genome sequences in FASTA (all segments)

    simple_filter_illumina_INFL.py ${bam} ${primers} ${params.max_depth} ${params.min_mapq} ${params.length} ${ref_genome_fasta}
    """
}
