process nextclade {
    tag "nextclade:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "nextstrain_lineage.csv"
    containerOptions "--volume ${params.nextclade_db_absolute_path_on_host}:/home/external_databases/nextclade_db"

    input:
    tuple val(sampleId), path(consensus_masked_sv_fa)

    output:
    tuple val(sampleId), path('nextstrain_lineage.json')
    tuple val(sampleId), path('nextclade_lineages/nextclade.cds_translation.S.fasta')

    script:
    """
    nextclade run --input-dataset /home/external_databases/nextclade_db/sars-cov-2.zip \
                  --output-csv nextstrain_lineage.csv \
                  --output-all nextclade_lineages \
                  ${consensus_masked_sv_fa}
    parse_nextclade_output_csv2json.py nextstrain_lineage.csv nextclade_lineages/nextclade.auspice.json nextstrain_lineage_sars.json full_genome

    # Combine two jsons in one
    jq -s "." nextstrain_lineage_*.json > nextstrain_lineage.json
    """
}
