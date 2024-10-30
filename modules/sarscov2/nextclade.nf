process nextclade {
    tag "nextclade:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "nextstrain_lineage.csv"
    containerOptions "--volume ${params.external_databases_path}:/home/external_databases/"

    input:
    tuple val(sampleId), path('consensus_masked_SV.fa'), path(ref_genome_with_index), val(QC_status)

    output:
    tuple val(sampleId), path("nextclade_lineages"), path('nextstrain_lineage.csv'), emit: to_pubdir
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
       touch nextstrain_lineage.csv
       mkdir nextclade_lineages
       touch nextclade_lineages/nextclade.cds_translation.S.fasta
       QC_status_exit="nie"
    else
      # This module need to handel both RSV and SARS, Influenza has a separate module  due to specific requirements
      # WE guess the correct nextclade file by checking fasta header
      HEADER=`head -1 ${ref_genome_with_index[final_index]} | tr -d ">"`
      if [[ "\${HEADER}" == "MN"* ]]; then
        NEXCLADE_FILE="sars-cov-2.zip"
      elif [[ \${HEADER} == "hRSV/A"* ]]; then
        NEXCLADE_FILE="RSV-A.zip"
      elif [[ \${HEADER} == "hRSV/B/"* ]]; then
        NEXCLADE_FILE="RSV-B.zip"
      fi

      nextclade run --input-dataset /home/external_databases/nextclade/\${NEXCLADE_FILE} \
                    --output-csv nextstrain_lineage.csv \
                    --output-all nextclade_lineages \
                    consensus_masked_SV.fa
      
      if [ -e nextclade_lineages/nextclade.cds_translation.S.fasta ]; then
        QC_status_exit="tak"
      else
        # tworzymy pusty plik ale modeller dostanie QC zeby sie nie wykonywal
        touch nextclade_lineages/nextclade.cds_translation.S.fasta
        QC_status_exit="nie"
      fi
    fi
    """
}
