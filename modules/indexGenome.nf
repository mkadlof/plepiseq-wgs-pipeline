process indexGenome {
    input:
    path reference_fasta

    output:
    tuple path("${reference_fasta}"), path("${reference_fasta}.fai")

    script:
    """
    samtools faidx ${reference_fasta}
    """
}
