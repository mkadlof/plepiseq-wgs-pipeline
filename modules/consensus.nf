process consensus {
    tag "Calculating consensus sequence for sample:\t$sampleId"
    // publishDir "${params.results_dir}/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(freebayes_fa), path(lofreq_fa), path(varscan_fa)

    output:
    tuple val(sampleId), path("consensus.fa")

    script:
    """
    cat ${freebayes_fa} ${lofreq_fa} ${varscan_fa} > tmp_to_consensus.fa
    mafft --auto --inputorder --quiet tmp_to_consensus.fa > tmp_to_consensus_aln.fa
    generate_consensus_sequence.py tmp_to_consensus_aln.fa
    mv output_consensus.fa consensus.fa
    """
}
