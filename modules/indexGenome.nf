process indexGenome {
    input:
    val(x)
    // path reference_fasta

    output:
    tuple path("genome.fasta"), path("genome.fasta.fai")

    script:
    """
    cp $x genome.fasta  
    samtools faidx genome.fasta
    """
}
