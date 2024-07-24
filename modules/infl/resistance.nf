process resistance {
    tag "resistance:${sampleId}"

    input:
    tuple val(sampleId), env("REF_GENOME_ID_MINI")
    tuple val(sampleId), path("nextalign_gene_*.translation.fasta")

    output:
    tuple val(sampleId), path("Mutation_report_*.txt")

    script:
    """
    # === This workaround because nextflow replace glob patterns with numbers ===
    # It goal is to rename symlinks to their targets names.
    for link in \$(find . -maxdepth 1 -type l); do
        target=\$(readlink "\$link")
        base_target=\$(basename "\$target")
        mv "\$link" "\$base_target"
    done
    # === End of workaround ===

    cat nextalign_gene_{NA,PA}.translation.fasta >> to_resistance.fasta
    TO_RESISTANCE='H1N1 H3N2 H5N1 H7N9 Yamagata Victoria'
    if [[ \${TO_RESISTANCE[@]} =~ \${REF_GENOME_ID_MINI} ]]; then
        analyze_infl_mutations.py to_resistance.fasta \${REF_GENOME_ID_MINI} analiza_opornosci /home/data/infl .
    fi
    """
}
