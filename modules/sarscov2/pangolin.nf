process pangolin {
    tag "pangolin:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "pangolin_lineage.csv"
    containerOptions "--volume ${params.pangolin_db_absolute_path_on_host}:/home/external_databases/pangolin"

    input:
    tuple val(sampleId), path(consensus_masked_sv_fa)

    output:
    tuple val(sampleId), path('pangolin.json')

    script:
    """
    pangolin --outfile pangolin_lineage.csv \
             --threads ${params.threads} \
             ${consensus_masked_sv_fa}
    parse_pangolin_output_csv2json.py pangolin_lineage.csv pangolin.json
    """
}
