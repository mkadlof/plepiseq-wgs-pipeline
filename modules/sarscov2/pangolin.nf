process pangolin {
    tag "pangolin:${sampleId}"
    container = params.main_image
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "pangolin_lineage.csv"
    containerOptions "--volume ${params.external_databases_path}:/home/external_databases/"

    input:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path(ref_genome_with_index), val(QC_status)

    output:
    tuple val(sampleId), path('pangolin.json'), emit:json
    script:
    """
    # Nie wiem jak ale pangolin po zamontowaniu external databases wie ze ma pobranego pangolin data
    if [[ ${QC_status} == "nie" || ${params.species} != "SARS-CoV-2" ]]; then
        touch pangolin_lineage.csv
        touch pangolin.json # wrong species or failed QC
        #TODO: add error message with upstream error
    else
        pangolin --outfile pangolin_lineage.csv \
               --threads ${params.threads} \
               output_consensus_masked_SV.fa
        parse_pangolin_output_csv2json.py pangolin_lineage.csv pangolin_sars.json
        jq -s "." pangolin_sars.json > pangolin.json
    fi
    """
}
