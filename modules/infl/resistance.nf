process resistance {
    tag "resistance:${sampleId}"
    container = params.main_image
    input:
    tuple val(sampleId), path("to_resistance.fasta"), path("sample_M2.fasta"), val(QC_status), val(SAMPLE_SUBTYPE)

    output:
    tuple val(sampleId), path("Mutation_report_*.txt"), emit: to_pubdir
    tuple val(sampleId), path("drug_resistance.json"), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch Mutation_report_dummy.txt
      touch drug_resistance.json
    else
  
      TO_RESISTANCE='H1N1 H3N2 H5N1 H7N9 Yamagata Victoria'
      if [[ \${TO_RESISTANCE[@]} =~ ${SAMPLE_SUBTYPE} ]]; then
          analyze_infl_mutations.py to_resistance.fasta ${SAMPLE_SUBTYPE} analiza_opornosci /home/data/infl .
          touch drug_resistance.json
      else
         # json with error message that for the provided subtype we do not determine drug resistance
         touch Mutation_report_dummy.txt
         touch drug_resistance.json
      fi
    fi
    """
}
