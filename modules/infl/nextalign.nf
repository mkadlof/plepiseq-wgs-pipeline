process nextalign {
    tag "nextalign:${sampleId}"
    container = params.main_image
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "*.fasta"

    input:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path(ref_genome_with_index), val(QC_status), val(SAMPLE_SUBTYPE)

    output:
    tuple val(sampleId), path("to_resistance.fasta"), path("sample_M2.fasta"), val(QC_status), val(SAMPLE_SUBTYPE)

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch to_resistance.fasta
      touch sample_M2.fasta
    else
      cat output_consensus_masked_SV.fa | awk '{if (substr(\$0, 1, 1)==">") { new_name=\$0; gsub("\\\\.", "_", new_name); gsub("/", "_", new_name);  gsub("_SV", "", new_name);  filename=("sample_"substr(new_name,2) ".fasta"); print \$0 >> filename } else {print toupper(\$0)  >> filename}}'
      
      NEXTALIGN_DB="/home/data/infl/nextalign"
      
      # For these subtypes we prepared custom files to get sequences of ALL proteins 
      TO_NEXTALIGN='H3N2 H5N1 H5N5 H5N6 H5N8 H7N9'
      if [[ \${TO_NEXTALIGN[@]} =~ ${SAMPLE_SUBTYPE} ]]; then
        for FILE in `ls sample_*.fasta`
        do
            SEGMENT=`basename \${FILE} ".fasta" | cut -d "_" -f3 `
            nextalign run -r \${NEXTALIGN_DB}/${SAMPLE_SUBTYPE}/\${SEGMENT}.fasta \
                          -m \${NEXTALIGN_DB}/${SAMPLE_SUBTYPE}/${SAMPLE_SUBTYPE}.gff \
                          -g \${SEGMENT} \
                          -O . \${FILE}

            if [[ \${SEGMENT} == 'NA' || \${SEGMENT} == 'PA' ]]; then
                cat nextalign_gene_\${SEGMENT}.translation.fasta >> to_resistance.fasta
            fi

            if [ \${SEGMENT} == 'MP' ]; then
                # Te jest glupie ale moj skrypt ma hard coded ta nazwe 
                cp sample_chr7_MP.fasta consensus_MP.fasta
                # I hard kodowny naglowek
                sed -i s"/chr7_MP_SV/MP/"g consensus_MP.fasta
                prep_M2.py ${SAMPLE_SUBTYPE} sample_M2.fasta \${NEXTALIGN_DB} . .
                nextalign run -r M2.fasta \
                              -m \${NEXTALIGN_DB}/${SAMPLE_SUBTYPE}/${SAMPLE_SUBTYPE}.gff \
                              -g M2 \
                              -O . consensus_M2.fasta
               cp nextalign_gene_M2.translation.fasta sample_M2.fasta
            fi
        done
      fi # koniec if-a na dostepnosc bazy nextalign
    fi
    """
}
