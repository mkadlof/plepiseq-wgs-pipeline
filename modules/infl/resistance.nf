process resistance {
    tag "resistance:${sampleId}"
    container = params.main_image
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path("to_resistance.fasta"), path("sample_M2.fasta"), val(QC_status), val(SAMPLE_SUBTYPE)

    output:
    tuple val(sampleId), path("drug_resistance.json"), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      STATUS="nie"
      if [ "${params.lan}" == "pl" ]; then
        ERR_MSG="Ten moduł został uruchomiony na próbce, która nie przeszła kontroli jakości."
      else
        ERR_MSG="This sample failed a QC analysis during an earlier phase of the analysis."
      fi
      analyze_infl_mutations.py --status \${STATUS} \
                                --error "\${ERR_MSG}" \
                                --output_json drug_resistance.json \
                                --lan ${params.lan}

    else
  
      TO_RESISTANCE='H1N1 H3N2 H5N1 H7N9 Yamagata Victoria'
      if [[ \${TO_RESISTANCE[@]} =~ ${SAMPLE_SUBTYPE} ]]; then
          analyze_infl_mutations.py --input_fasta to_resistance.fasta \
                                    --subtype ${SAMPLE_SUBTYPE} \
                                    --sample_name ${sampleId} \
                                    --data_path /home/data/infl \
                                    --output_path . \
                                    --status "tak" \
                                    --output_json drug_resistance.json \
                                    --lan ${params.lan}
      else
         STATUS="nie"
         if [ "${params.lan}" == "pl" ]; then
           ERR_MSG="Ten moduł działa tylko dla następujących podtypów wirusa grypy: H1N1, H3N2, H5N1, H7N9, Yamagata, Victoria. Natomiast podtyp próbki to: ${SAMPLE_SUBTYPE}"
         else
           ERR_MSG="This module works only for following subtypes: H1N1 H3N2 H5N1 H7N9 Yamagata Victoria. Sample, however is: ${SAMPLE_SUBTYPE}"
         fi
         analyze_infl_mutations.py --status \${STATUS} \
                                    --error "\${ERR_MSG}" \
                                    --output_json drug_resistance.json \
                                    --lan ${params.lan}
      fi
    fi
    """
}
