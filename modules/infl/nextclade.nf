process nextclade {
    
    tag "nextclade for ${SEGMENT}:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "nextstrain_lineage.csv"
    
    container  = params.main_image
    containerOptions "--volume ${params.external_databases_path}:/home/external_databases/"

    input:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path(ref_genome_with_index), val(QC_status), val(SAMPLE_SUBTYPE)

    output:
    tuple val(sampleId), path("nextclade_lineage_*.csv"), emit: to_pubdir
    tuple val(sampleId), path("nextclade.json"), emit: json

    script:
    """
    # we split final genome into segments
    if [ ${QC_status} == "nie" ]; then
      touch nextclade_lineage_dummy.csv
      touch nextclade.json
    else

      cat output_consensus_masked_SV.fa | awk '{if (substr(\$0, 1, 1)==">") { new_name=\$0; gsub("\\\\.", "_", new_name); gsub("/", "_", new_name);  gsub("_SV", "", new_name);  filename=("sample_"substr(new_name,2) ".fasta"); print \$0 >> filename } else {print toupper(\$0)  >> filename}}'
   
      # For this subtypes nextclade provides a database
      KNOWN='H1N1 H3N2 Yamagata Victoria'
      if [[ \${KNOWN[@]} =~ ${SAMPLE_SUBTYPE} ]]; then
        nextclade run \
            --input-dataset /home/external_databases/nextclade/${SAMPLE_SUBTYPE}_HA.zip \
            --output-csv nextclade_lineage_HA.csv \
            --output-all nextclade_linages_HA \
            sample_chr4_HA.fasta
        # for all this subtypes, save, Yamagata, we can also analyze NA
        if [ ${SAMPLE_SUBTYPE} != 'Yamagata' ]; then
            nextclade run --input-dataset /home/external_databases/nextclade/${SAMPLE_SUBTYPE}_NA.zip \
                --output-csv nextclade_lineage_NA.csv \
                --output-all nextclade_linages_NA \
                sample_chr6_NA.fasta
      
        fi
        # placeholder for json 
        touch nextclade.json
      else
        # json with message that for provided subtype there is no nextclade database
        touch nextclade_lineage_dummy.csv
        touch nextclade.json    
      fi # koniec if-a na subtypes w nextclade
    fi # koniec if-a na QC
    """
}
