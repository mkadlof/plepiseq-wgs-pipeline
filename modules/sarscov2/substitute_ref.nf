process substitute_ref_genome {
    // This process is used to substitute reference genome from second stept of genome generation
    // in nanopore to the originally used reference genome and to integrate nanopore and illumina I/O
    tag "genome_substitution:${sampleId}"
    container  = params.main_image
    input:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path(ref_genome_with_index), val(QC_status), path('original.fasta')
    output:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path('original_reference.*'), val(QC_status), emit: fasta_refgenome_and_qc
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), val(QC_status), emit: fasta_and_qc // for novel snpeff
    script:
    """
   if [ ${QC_status} == "nie" ]; then
     touch output_consensus_masked_SV.fa
     touch original_reference.fasta
   else
     cp original.fasta original_reference.fasta
     bwa index original_reference.fasta
   fi
   """
}
