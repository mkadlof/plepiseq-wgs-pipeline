process coinfection_analysis_illumina {
    tag "coinfection_analysis:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "*allele_usage_histogram.txt"

    input:
    tuple val(sampleId), path(detected_variants_varscan_coinfection_txt), val(QC_status)

    output:
    tuple val(sampleId), path("${sampleId}_allele_usage_histogram.txt"), emit: to_pubdir
    tuple val(sampleId), path('custom_coinfection_analysis.json'), emit: json

    script:
    """
    # The coinfection analysis is based on similarity to 3 known samples from EQA2023.
    # The new version of the script allows specifying any number of samples that are known to be
    # coinfected. Simply add additional files as positional arguments.

    if [ ${QC_status} == "nie" ]; then
      touch ${sampleId}_allele_usage_histogram.txt
      coinfection_status="nie"
      echo -e "{\\"coinfection_status\\":\\"\${coinfection_status}\\"}" >> custom_coinfection_analysis.json
    else
      RESULTS=(`predict_coinfection_illumina.py ${detected_variants_varscan_coinfection_txt} \
                                      ${sampleId} \
                                      /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS2.09_coinfections.txt \
                                      /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS2.17_coinfections.txt \
                                      /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS2.32_coinfections.txt`)
      cp allele_usage_histogram.txt ${sampleId}_allele_usage_histogram.txt
      coinfection_status="tak"
      coinfection_result=\${RESULTS[0]}
      coinfection_pvalue=`echo \${RESULTS[1]} | awk '{printf "%.2f", \$0}'`
      coinfection_histogram_file="${params.results_dir}/${sampleId}/${sampleId}_allele_usage_histogram.txt"

      echo -e "{\\"coinfection_status\\":\\"\${coinfection_status}\\",
                \\"coinfection_result\\":\\"\${coinfection_result}\\",
                \\"coinfection_pvalue\\":\${coinfection_pvalue},
                \\"coinfection_histogram_file\\":\\"\${coinfection_histogram_file}\\"}" >> custom_coinfection_analysis.json

    fi

    """
}

process coinfection_analysis_nanopore {
    // The only difference with respect to illumina is that we use different "reference" files
    tag "coinfection_analysis:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "*allele_usage_histogram.txt"

    input:
    tuple val(sampleId), path(detected_variants_varscan_coinfection_txt), val(QC_status)

    output:
    tuple val(sampleId), path("${sampleId}_allele_usage_histogram.txt"), emit: to_pubdir
    tuple val(sampleId), path('custom_coinfection_analysis.json'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch ${sampleId}_allele_usage_histogram.txt
      coinfection_status="nie"
      echo -e "{\\"coinfection_status\\":\\"\${coinfection_status}\\"}" >> custom_coinfection_analysis.json
    else
      RESULTS=(`predict_coinfection_illumina.py ${detected_variants_varscan_coinfection_txt} \
                                      ${sampleId} \
                                      /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS1.04_contamination.txt \
                                      /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS1.05_contamination.txt`)
      cp allele_usage_histogram.txt ${sampleId}_allele_usage_histogram.txt

      coinfection_status="tak"
      coinfection_result=\${RESULTS[0]}
      coinfection_pvalue=`echo \${RESULTS[1]} | awk '{printf "%.2f", \$0}'`
      coinfection_histogram_file="${params.results_dir}/${sampleId}/${sampleId}_allele_usage_histogram.txt"

      echo -e "{\\"coinfection_status\\":\\"\${coinfection_status}\\",
                \\"coinfection_result\\":\\"\${coinfection_result}\\",
                \\"coinfection_pvalue\\":\${coinfection_pvalue},
                \\"coinfection_histogram_file\\":\\"\${coinfection_histogram_file}\\"}" >> custom_coinfection_analysis.json
    fi
    """
} 
