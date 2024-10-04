process nextclade {
    tag "nextclade:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "nextstrain_lineage.csv"
    containerOptions "--volume ${params.nextclade_db_absolute_path_on_host}:/home/external_databases/nextclade_db"

    input:
    tuple val(sampleId), env("REF_GENOME_ID_MINI"), path("consensus_*.fasta")

    output:
    tuple val(sampleId), path("nextstrain_lineage_HA.csv"), path("nextclade_lineages_HA"), path('nextstrain_lineage_NA.csv'), path("nextclade_lineages_NA")

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

    touch nextstrain_lineage_HA.csv
    touch nextclade_lineages_HA
    touch nextstrain_lineage_NA.csv
    touch nextclade_lineages_NA

    echo "REF_GENOME_ID_MINI: \${REF_GENOME_ID_MINI}"

    KNOWN='H1N1 H3N2 Yamagata Victoria'
    if [[ \${KNOWN[@]} =~ \${REF_GENOME_ID_MINI} ]]; then
        nextclade run \
            --input-dataset /home/external_databases/nextclade_db/\${REF_GENOME_ID_MINI}_HA.zip \
            --output-csv nextstrain_lineage_HA.csv \
            --output-all nextclade_linages_HA \
            consensus_HA.fasta

        if [ \${REF_GENOME_ID_MINI} != 'Yamagata' ]; then
            nextclade run --input-dataset /home/external_databases/nextclade_db/\${REF_GENOME_ID_MINI}_NA.zip \
                --output-csv nextstrain_lineage_NA.csv \
                --output-all nextclade_linages_NA \
                consensus_NA.fasta
        fi
    fi
    """
}
