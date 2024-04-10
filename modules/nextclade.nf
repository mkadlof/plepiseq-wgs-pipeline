process nextclade {
    publishDir "results/${sampleId}", mode: 'symlink'
    containerOptions "--volume ${params.nextclade_db_absolute_path_on_host}:/home/external_databases/nextclade_db"

    input:
    tuple val(sampleId), path(consensus_masked_sv_fa)

    output:
    tuple val(sampleId), path("nextclade_lineages"), path('nextstrain_lineage.csv')
    tuple val(sampleId), path('nextclade_lineages/nextclade.cds_translation.S.fasta')

    script:
    """
    nextclade run --input-dataset /home/external_databases/nextclade_db/sars-cov-2.zip \
                  --output-csv nextstrain_lineage.csv \
                  --output-all nextclade_lineages \
                  ${consensus_masked_sv_fa}
    """
}