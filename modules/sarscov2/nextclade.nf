process nextclade {
    tag "nextclade:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "nextstrain_lineage.csv"
    containerOptions "--volume ${params.external_databases_path}:/home/external_databases/"

    input:
    tuple val(sampleId), path('consensus_masked_SV.fa'), path(ref_genome_with_index), val(QC_status)

    output:
    tuple val(sampleId), path('nextstrain_lineage.json'), emit: json
    tuple val(sampleId), path('nextclade_lineages/nextclade.cds_translation.S.fasta'), env(QC_status_exit), emit:to_modeller

    script:
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }
    """
    if [ ${QC_status} == "nie" ]; then
       mkdir nextclade_lineages
       touch nextclade_lineages/nextclade.cds_translation.S.fasta
       QC_status_exit="nie"
      cat <<EOF > nextstrain_lineage.json
[
    {
      "status": "nie",
      "database_name": "Nextclade",
      "error_message": "QC failed: an error occurred in a prior processing step."
    }
]
EOF
    else
      # This module need to handel both RSV and SARS, Influenza has a separate module  due to specific requirements
      # WE guess the correct nextclade file by checking fasta header
      HEADER=`head -1 ${ref_genome_with_index[final_index]} | tr -d ">"`
      if [[ "\${HEADER}" == "MN"* ]]; then
        NEXCLADE_FILE="sars-cov-2.zip"
      elif [[ \${HEADER} == "hRSV/A"* ]]; then
        NEXCLADE_FILE="RSV_A.zip"
      elif [[ \${HEADER} == "hRSV/B/"* ]]; then
        NEXCLADE_FILE="RSV_B.zip"
    fi
    nextclade run --input-dataset /home/external_databases/nextclade/\${NEXCLADE_FILE} \
                    --output-csv nextstrain_lineage.csv \
                    --output-all nextclade_lineages \
                    consensus_masked_SV.fa

     if [ -e nextclade_lineages/nextclade.cds_translation.S.fasta ]; then
        parse_nextclade_output_csv2json.py nextstrain_lineage.csv nextclade_lineages/nextclade.auspice.json nextstrain_lineage_sars.json full_genome

        # Combine two jsons in one
        jq -s "." nextstrain_lineage_*.json > nextstrain_lineage.json
        QC_status_exit="tak"
      else
        # tworzymy pusty plik ale modeller dostanie QC zeby sie nie wykonywal
        touch nextclade_lineages/nextclade.cds_translation.S.fasta
        QC_status_exit="nie"
        touch nextstrain_lineage.json
        mkdir nextclade_lineages || true
        touch nextclade_lineages/nextclade.cds_translation.S.fasta
        QC_status_exit="nie"
        #TODO: add error_mesage to json output file.
      fi
    fi
    """
}
