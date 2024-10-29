process consensus {
    tag "consensus:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "consensus_*.fasta"

    input:
    tuple val(sampleId), path(masked_ref_genome_fa), path(varscan_fa), path(freebayes_fa), val(QC_status), path(lofreq_fa)

    output:
    tuple val(sampleId), path("consensus.fasta"), val(QC_status), emit: single_fasta
    tuple val(sampleId), path("consensus_*.fasta"), val(QC_status), emit: multiple_fastas
    tuple val(sampleId), path("consensus.json"), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch consensus.fasta
      touch consensus_A.fasta
      touch consensus.json

    else

      make_consensus.py ${masked_ref_genome_fa} ${freebayes_fa} ${lofreq_fa} ${varscan_fa}
      # get the total length and number of Ns in the consensus

      TOTAL_LENGTH=\$(grep -v '>' consensus.fasta | wc -c)
      NUMBER_OF_N=\$(grep -v '>' consensus.fasta | grep N -o | wc -l)

      echo -e '{"total_length_value": \${TOTAL_LENGTH}\n"number_of_Ns_value": \${NUMBER_OF_N}}' >> consensus.json 
    fi
    """
}
