process indexGenome {

    input:
    path ref_genome

    output:
    path("${ref_genome}.fai")

    script:
    """
    samtools faidx ${ref_genome}
    """
}