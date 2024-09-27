process consensus {
    tag "consensus:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "consensus_*.fasta"

    input:
    tuple val(sampleId), path(ref_genome_fa), path(freebayes_fa), path(lofreq_fa), path(varscan_fa)

    output:
    tuple val(sampleId), path("consensus.fasta")
    tuple val(sampleId), path("consensus_*.fasta")

    script:
    """
    make_consensus.py ${ref_genome_fa} ${freebayes_fa} ${lofreq_fa} ${varscan_fa}
    """
}
