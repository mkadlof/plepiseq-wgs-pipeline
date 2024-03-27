process nextclade {
    publishDir "results/${sampleId}", mode: 'symlink'
    containerOptions "--volume ${params.nextclade_db_absolute_path_on_host}:/home/SARS-CoV2/nextclade_db"

    input:
    tuple val(sampleId), path(consensus_masked_sv_fa)

    output:
    tuple val(sampleId), path("nextclade_lineages"), path('nextstrain_lineage.csv')

    script:
    """
    nextclade run --input-dataset /home/SARS-CoV2/nextclade_db/sars-cov-2.zip \
                  --output-csv nextstrain_lineage.csv \
                  --output-all nextclade_lineages \
                  ${consensus_masked_sv_fa}
    """
}