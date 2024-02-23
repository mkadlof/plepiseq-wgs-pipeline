process consensus {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(freebayes_fa), path(lofreq_fa), path(varscan_fa)

    output:
    tuple val(sampleId), path("consensus.fa")

    script:
    """
    cat ${freebayes_fa} ${lofreq_fa} ${varscan_fa} > tmp_to_consensus.fa
    mafft --auto --inputorder --quiet tmp_to_consensus.fa > tmp_to_consensus_aln.fa
    gen_consensus_seq_final.py tmp_to_consensus_aln.fa
    mv output_consensus.fa consensus.fa
    """
}