process consensus_illumina {
    // For illumina this module provides sample sequnce PRIOR to introduction of SVs
    tag "consensus:${sampleId}"
    container  = params.main_image
    // publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "consensus_*.fasta"

    input:
    tuple val(sampleId), path(masked_ref_genome_fa), path(varscan_fa), path(freebayes_fa), val(QC_status), path(lofreq_fa)

    output:
    // tuple val(sampleId), path("consensus.fasta"), val(QC_status), emit: single_fasta
    tuple val(sampleId), path("consensus_*.fasta"), val(QC_status), emit: multiple_fastas
    // tuple val(sampleId), path("consensus.json"), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch consensus.fasta
      touch consensus_dummy_segment.fasta
      ERR_MSG="Failed QC"
      parse_make_consensus.py --status "nie" --error "\${ERR_MSG}" -o consensus.json
    else
      make_consensus.py ${masked_ref_genome_fa} ${freebayes_fa} ${lofreq_fa} ${varscan_fa}
      # prepare json for this step
      ls consensus_*.fasta | tr " " "\\n" >> list_of_fasta.txt
      parse_make_consensus.py --status "tak" -o consensus.json --input_fastas list_of_fasta.txt --output_path "${params.results_dir}/${sampleId}"

    fi
    """
}

process consensus_nanopore {
    // For nanopore this module provides FINAL sequence for a sample
    tag "consensus:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "consensus_*.fasta"

    input:
    tuple val(sampleId), path(masked_ref_genome_fa), path(sample_genome), val(QC_status), path('genome.fasta')

    output:
    // tuple val(sampleId), path("consensus.fasta"), val(QC_status), emit: single_fasta
    tuple val(sampleId), path("consensus_*.fasta"), val(QC_status), emit: multiple_fastas
    tuple val(sampleId), path('consensus.fasta'), path('ref_genome.*'), val(QC_status), emit: fasta_refgenome_and_qc
    tuple val(sampleId), path("consensus.json"), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch consensus.fasta
      touch consensus_dummy_segment.fasta
      touch ref_genome.fasta
      touch ref_genome.fasta.fai
      ERR_MSG="Failed QC"
      parse_make_consensus.py --status "nie" --error "\${ERR_MSG}" -o consensus.json
    else
      make_consensus_nanopore.py ${masked_ref_genome_fa} ${sample_genome}

      # This step is required only for integration with downstream illumina modules
      cp genome.fasta ref_genome.fasta
      bwa index ref_genome.fasta
      # prepare json for this step
      ls consensus_*.fasta | tr " " "\\n" >> list_of_fasta.txt
      parse_make_consensus.py --status "tak" -o consensus.json --input_fastas list_of_fasta.txt --output_path "${params.results_dir}/${sampleId}"
      sed -i s'|_SV"|"|'g consensus.json
    
    fi
    """
}
