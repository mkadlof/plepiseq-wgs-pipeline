process consensusMasking {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(lowcoverage_masked_fa), path(consensus_fa)

    output:
    tuple val(sampleId), path("consensus_masked.fa")

    script:
    """
    cat lowcoverage_masked.fa consensus.fa >tmp_consensus.fa
    mafft --auto --inputorder --quiet tmp_consensus.fa > tmp_consensus_aln.fa
    get_N.py tmp_consensus_aln.fa
    mv output_consensus_masked.fa consensus_masked.fa
    """
}