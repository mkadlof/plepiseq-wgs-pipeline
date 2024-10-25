process nextclade {
    tag "nextclade:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "nextstrain_lineage.csv"
    containerOptions "--volume ${params.nextclade_db_absolute_path_on_host}:/home/external_databases/nextclade_db"

    input:
    tuple val(sampleId), path('consensus_masked_SV.fa'), path(ref_genome_with_index), val(QC_status)

    output:
    tuple val(sampleId), path("nextclade_lineages"), path('nextstrain_lineage.csv'), emit: to_pubdir
    tuple val(sampleId), path('nextclade_lineages/nextclade.cds_translation.S.fasta'), optional: true, emit:to_modeller

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
       touch nextstrain_lineage.csv
    else
      # We assueme that we dont't now which pipeline uses his module and guess that by checking reference
      # sequence header
      HEADER=`head -1 ${ref_genome_with_index[final_index]} | tr -d ">"`
      if [[ "\${HEADER}" == "MN"* ]]; then
        NEXCLADE_FILE="sars-cov-2.zip"
      elif [[ \${HEADER} == "hRSV/A/"* ]]; then
        NEXCLADE_FILE="RSV-A.zip"
      elif [[ \${HEADER} == "hRSV/B/"* ]]; then
        NEXCLADE_FILE="RSV-B.zip"
      fi

      nextclade run --input-dataset /home/external_databases/nextclade_db/\${NEXCLADE_FILE} \
                    --output-csv nextstrain_lineage.csv \
                    --output-all nextclade_lineages \
                    consensus_masked_SV.fa
    fi
    """
}
