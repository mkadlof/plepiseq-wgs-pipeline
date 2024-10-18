process consensus {
    tag "consensus:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "consensus_*.fasta"

    input:
    tuple val(sampleId), path(ref_genome_fa), path(freebayes_fa), path(lofreq_fa), path(varscan_fa)

    output:
    tuple val(sampleId), path("consensus.fasta")
    tuple val(sampleId), path("consensus_*.fasta")
    tuple val(sampleId), path("consensus.json")

    script:
    """
    make_consensus.py ${ref_genome_fa} ${freebayes_fa} ${lofreq_fa} ${varscan_fa}

    # get the total length and number of Ns in the consensus

    TOTAL_LENGTH=\$(grep -v '>' consensus.fasta | wc -c)
    NUMBER_OF_N=\$(grep -v '>' consensus.fasta | grep N -o | wc -l)

    cat > consensus.json <<EOF
    {
        "total_length_value": \${TOTAL_LENGTH},
        "number_of_Ns_value": \${NUMBER_OF_N}
    }
    EOF
    """
}
