process nextclade {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(consensus_masked_sv_fa)

    output:
    tuple val(sampleId), path("nextclade_lineages"), path('nextstrain_lineage.csv')

    script:
    """
    nextclade run --input-dataset ${params.nextclade_db_dir} \
                  --output-csv nextstrain_lineage.csv \
                  --output-all nextclade_lineages \
                  ${consensus_masked_sv_fa}
    """
}