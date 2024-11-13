process consensus_illumina {
    tag "consensus:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "consensus_*.fasta"

    input:
    tuple val(sampleId), path(masked_ref_genome_fa), path(varscan_fa), path(freebayes_fa), val(QC_status), path(lofreq_fa)

    output:
    tuple val(sampleId), path("consensus.fasta"), val(QC_status), emit: single_fasta
    tuple val(sampleId), path("consensus_*.fasta"), val(QC_status), emit: multiple_fastas
    tuple val(sampleId), path("consensus.json"), path("list_of_fasta_files.txt"), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch consensus.fasta
      touch consensus_dummy_segment.fasta
      touch consensus.json
      ls consensus_*.fasta > list_of_fasta_files.txt
    else
      make_consensus.py ${masked_ref_genome_fa} ${freebayes_fa} ${lofreq_fa} ${varscan_fa}
      # get the total length and number of Ns in the consensus

      TOTAL_LENGTH=\$(grep -v '>' consensus.fasta | wc -c)
      NUMBER_OF_N=\$(grep -v '>' consensus.fasta | grep N -o | wc -l)

      echo -e "{\\"total_length_value\\": \${TOTAL_LENGTH},\n\\"number_of_Ns_value\\": \${NUMBER_OF_N}}" >> consensus.json

      ls consensus_*.fasta > list_of_fasta_files.txt
    fi
    """
}

process consensus_nanopore {
    tag "consensus:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "consensus_*.fasta"

    input:
    tuple val(sampleId), path(masked_ref_genome_fa), path(sample_genome), val(QC_status)

    output:
    stdout
    // tuple val(sampleId), path("consensus.fasta"), val(QC_status), emit: single_fasta
    // tuple val(sampleId), path("consensus_*.fasta"), val(QC_status), emit: multiple_fastas
    // tuple val(sampleId), path("consensus.json"), path("list_of_fasta_files.txt"), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch consensus.fasta
      touch consensus_dummy_segment.fasta
      touch consensus.json
      ls consensus_*.fasta > list_of_fasta_files.txt
    else
      make_consensus_nanopore.py ${masked_ref_genome_fa} ${sample_genome}
      # get the total length and number of Ns in the consensus

      #TOTAL_LENGTH=\$(grep -v '>' consensus.fasta | wc -c)
      #NUMBER_OF_N=\$(grep -v '>' consensus.fasta | grep N -o | wc -l)

      # echo -e "{\\"total_length_value\\": \${TOTAL_LENGTH},
      # \\"number_of_Ns_value\\": \${NUMBER_OF_N}}" >> consensus.json

      #ls consensus_*.fasta > list_of_fasta_files.txt
    fi
    """
}
