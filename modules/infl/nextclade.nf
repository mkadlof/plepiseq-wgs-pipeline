process nextclade {
    tag "nextclade:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "nextstrain_lineage.csv"
    containerOptions "--volume ${params.nextclade_db_absolute_path_on_host}:/home/external_databases/nextclade_db"

    input:
    tuple val(sampleId), env("REF_GENOME_ID_MINI"), path("consensus_*.fasta")

    output:
    tuple val(sampleId), path('nextstrain_lineage.json')
    tuple val(sampleId), path('nextstrain_lineage.json')

    script:
    """
    set -x

    # === This workaround because nextflow replace glob patterns with numbers ===
    # It goal is to rename symlinks to their targets names.
    for link in \$(find . -maxdepth 1 -type l); do
        target=\$(readlink "\$link")
        base_target=\$(basename "\$target")
        mv "\$link" "\$base_target"
    done
    # === End of workaround ===

    touch nextstrain_lineage_HA.csv
    mkdir nextclade_lineages_HA
    touch nextstrain_lineage_NA.csv
    mkdir nextclade_lineages_NA

    echo "REF_GENOME_ID_MINI: \${REF_GENOME_ID_MINI}"

    # Default json file
    generate_empty_json() {
        local segment="\$1"
        cat <<EOF > "nextstrain_lineage_\${segment}.json"
{
    "status": "nie",
    "database_name": "Nextclade",
    "error_message": "Nextclade was skip for \${REF_GENOME_ID_MINI} in sample ${sampleId} for segment \${segment}",
    "sequence_source": "\${segment}"
}
EOF
    }

    KNOWN='H1N1 H3N2 Yamagata Victoria'
    if [[ \${KNOWN[@]} =~ \${REF_GENOME_ID_MINI} ]]; then
        nextclade run \
            --input-dataset /home/external_databases/nextclade_db/\${REF_GENOME_ID_MINI}_HA.zip \
            --output-csv nextstrain_lineage_HA.csv \
            --output-all nextclade_lineages_HA \
            consensus_HA.fasta
            if [[ -s nextstrain_lineage_HA.csv ]]; then
                parse_nextclade_output_csv2json.py nextstrain_lineage_HA.csv nextclade_lineages_HA/nextclade.auspice.json nextstrain_lineage_HA.json HA
            fi
        if [ \${REF_GENOME_ID_MINI} != 'Yamagata' ]; then
            nextclade run --input-dataset /home/external_databases/nextclade_db/\${REF_GENOME_ID_MINI}_NA.zip \
                --output-csv nextstrain_lineage_NA.csv \
                --output-all nextclade_lineages_NA \
                consensus_NA.fasta
            if [[ -s nextstrain_lineage_NA.csv ]]; then
                parse_nextclade_output_csv2json.py nextstrain_lineage_NA.csv nextclade_lineages_NA/nextclade.auspice.json nextstrain_lineage_NA.json NA
            fi
        fi
    fi

    # If files where not created by Nextclade then create empty ones.
    if [[ ! -f "nextstrain_lineage_HA.json" ]]; then
        generate_empty_json "HA"
    fi

    if [[ ! -f "nextstrain_lineage_NA.json" ]]; then
        generate_empty_json "NA"
    fi

    # Combine two jsons in one
    jq -s "." nextstrain_lineage_HA.json nextstrain_lineage_NA.json > nextstrain_lineage.json
    """
}
