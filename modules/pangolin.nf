process pangolin {
    tag "Predicting PANGO lineage for for sample:\t$sampleId"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "pangolin_lineage.csv"
    containerOptions "--volume ${params.pangolin_db_absolute_path_on_host}:/home/external_databases/pangolin"

    input:
    tuple val(sampleId), path(consensus_masked_sv_fa)

    output:
    tuple val(sampleId), path('pangolin_lineage.csv')

    script:
    """
    pangolin --outfile pangolin_lineage.csv \
             --threads ${params.threads} \
             ${consensus_masked_sv_fa}
    """
}
