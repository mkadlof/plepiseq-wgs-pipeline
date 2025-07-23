process nextclade {
    tag "nextclade:${sampleId}"
    // publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "nextstrain_lineage.csv"
    container  = params.main_image
    containerOptions "--volume ${params.external_databases_path}:/home/external_databases/"
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path(ref_genome_with_index), val(QC_status), val(SAMPLE_SUBTYPE)

    output:
    tuple val(sampleId), path("nextclade_lineage.json"), emit: json
    tuple val(sampleId), path('dummy.fasta'), val(QC_status), emit:to_modeller

    script:
    """
    set -x
    touch dummy.fasta
    
 # Default json file
    generate_empty_json() {
    local segment="\$1"
    local lang="\$2"  # 'pl' or 'en'

    if [[ "\$lang" == "pl" ]]; then
        err_msg="Nextclade został pominięty dla podtypu ${SAMPLE_SUBTYPE} w próbce ${sampleId} dla segmentu \${segment}"
    else
        err_msg="Nextclade was skipped for ${SAMPLE_SUBTYPE} in sample ${sampleId} for segment \${segment}"
    fi

    cat <<EOF > "nextclade_lineage_\${segment}.json"
{
    "status": "nie",
    "database_name": "Nextclade",
    "error_message": "\${err_msg}",
    "sequence_source": "\${segment}"
}
EOF
}

    # we split final genome into segments
    if [ ${QC_status} == "nie" ]; then

        generate_empty_json "HA" ${params.lan}
        generate_empty_json "NA" ${params.lan}

        jq -s "." nextclade_lineage_HA.json nextclade_lineage_NA.json > nextclade_lineage.json
    else
        cat output_consensus_masked_SV.fa | awk '{if (substr(\$0, 1, 1)==">") { new_name=\$0; gsub("\\\\.", "_", new_name); gsub("/", "_", new_name);  gsub("_SV", "", new_name);  filename=("sample_"substr(new_name,2) ".fasta"); print \$0 >> filename } else {print toupper(\$0)  >> filename}}'
        KNOWN='H1N1 H3N2 Yamagata Victoria'
        if [[ \${KNOWN[@]} =~ ${SAMPLE_SUBTYPE} ]]; then
            nextclade run \
                --input-dataset /home/external_databases/nextclade/${SAMPLE_SUBTYPE}_HA.zip \
                --output-csv nextclade_lineage_HA.csv \
                --output-all nextclade_lineage_HA_dir \
                sample_chr4_HA.fasta
            if [[ -s nextclade_lineage_HA.csv ]]; then
                parse_nextclade_output_csv2json.py --input nextclade_lineage_HA.csv \
                                           --input2 nextclade_lineage_HA_dir/nextclade.auspice.json \
                                           --output nextclade_lineage_HA.json \
                                           --sequence_source HA \
                                           --lan ${params.lan}

            fi
            # for all this subtypes, save, Yamagata, we can also analyze NA
            if [ ${SAMPLE_SUBTYPE} != 'Yamagata' ]; then
                nextclade run --input-dataset /home/external_databases/nextclade/${SAMPLE_SUBTYPE}_NA.zip \
                    --output-csv nextclade_lineage_NA.csv \
                    --output-all nextclade_lineage_NA_dir \
                    sample_chr6_NA.fasta
                if [[ -s nextclade_lineage_NA.csv ]]; then
                    parse_nextclade_output_csv2json.py --input nextclade_lineage_NA.csv \
                                           --input2 nextclade_lineage_NA_dir/nextclade.auspice.json \
                                           --output nextclade_lineage_NA.json \
                                           --sequence_source NA \
                                           --lan ${params.lan}

                fi
            fi
        fi

        # If files where not created by Nextclade then create empty ones.
        if [[ ! -f "nextclade_lineage_HA.json" ]]; then
            generate_empty_json "HA" ${params.lan}
        fi

        if [[ ! -f "nextclade_lineage_NA.json" ]]; then
            generate_empty_json "NA" ${params.lan}
        fi

        # Combine two jsons in one
        jq -s "." nextclade_lineage_HA.json nextclade_lineage_NA.json > nextclade_lineage.json
    fi
    """
}
