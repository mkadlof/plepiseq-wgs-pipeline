process pangolin {
    tag "pangolin:${sampleId}"
    container = params.main_image
    cpus { params.threads > 3 ? 3 : params.threads }
    memory "20 GB"
    // publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "pangolin_lineage.csv"
    containerOptions "--volume ${params.external_databases_path}:/home/external_databases/"

    input:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path(ref_genome_with_index), val(QC_status)

    output:
    tuple val(sampleId), path('pangolin.json'), emit:json
    script:
    """
    # Nie wiem jak ale pangolin po zamontowaniu external databases wie ze ma pobranego pangolin data

    if [[ ${QC_status} == "nie" || ${params.species} != "SARS-CoV-2" ]]; then
        if [[ ${params.species} == "SARS-CoV-2" ]]; then
            if [ "${params.lan}" == "pl" ]; then
                ERR_MSG="Ten moduł został uruchomiony na próbce, która nie przeszła kontroli jakości."
            else
                ERR_MSG="This sample failed a QC analysis during an earlier phase of the analysis."
            fi
        elif [[ ${params.species} != "SARS-CoV-2" ]]; then
            if [ "${params.lan}" == "pl" ]; then
                ERR_MSG="Ten moduł jest przeznaczony wylacznie do analizy wirusa SARS-CoV-2."
            else
                ERR_MSG="For organisms other than SARS-CoV-2, the Pangolin database is not queried."
            fi
        fi
        touch pangolin_lineage.csv

        STATUS="nie"
        DATABASE_NAME="Pangolin"
        SEQUENCE_SOURCE="full_genome"
        echo -e "[{\\"status\\":\\"\${STATUS}\\",
                   \\"sequence_source\\":\\"\${SEQUENCE_SOURCE}\\",
                   \\"database_name\\":\\"\${DATABASE_NAME}\\",
                   \\"error_message\\":\\"\${ERR_MSG}\\"}]" >> pangolin.json

    else
        pangolin --outfile pangolin_lineage.csv \
                  --threads ${task.cpus} \
                  output_consensus_masked_SV.fa

        parse_pangolin_output_csv2json.py --input pangolin_lineage.csv \
                                          --output pangolin.json \
                                          --lan ${params.lan}
    fi
    """
}
