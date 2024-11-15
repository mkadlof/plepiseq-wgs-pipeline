process coinfection_analysis_illumina {
    tag "coinfection_analysis:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy'

    input:
    tuple val(sampleId), path(detected_variants_varscan_coinfection_txt), val(QC_status)

    output:
    tuple val(sampleId), path("${sampleId}_alternative_alleles_frequencies.png"), path("${sampleId}_coinfection_summary.txt"), emit: to_pubdir
    tuple val(sampleId), path('custom_coinfection_analysis.json'), emit: json

    script:
    """
    # The coinfection analysis is based on similarity to 3 known samples from EQA2023.
    # The new version of the script allows specifying any number of samples that are known to be
    # coinfected. Simply add additional files as positional arguments.

    touch custom_coinfection_analysis.json # Placeholder for json
    if [ ${QC_status} == "nie" ]; then
      touch ${sampleId}_coinfection_summary.txt
      touch ${sampleId}_alternative_alleles_frequencies.png
    else
      predict_coinfection_illumina.py ${detected_variants_varscan_coinfection_txt} \
                                      ${sampleId} \
                                      /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS2.09_coinfections.txt \
                                      /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS2.17_coinfections.txt \
                                      /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS2.32_coinfections.txt
    fi

    """
}

process coinfection_analysis_nanopore {
    // The only difference with respect to illumina is that we use different "reference" files
    tag "coinfection_analysis:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy'

    input:
    tuple val(sampleId), path(detected_variants_varscan_coinfection_txt), val(QC_status)

    output:
    tuple val(sampleId), path("${sampleId}_alternative_alleles_frequencies.png"), path("${sampleId}_coinfection_summary.txt"), emit: to_pubdir
    tuple val(sampleId), path('custom_coinfection_analysis.json'), emit: json

    script:
    """
    touch custom_coinfection_analysis.json # Placeholder for json
    if [ ${QC_status} == "nie" ]; then
      touch ${sampleId}_coinfection_summary.txt
      touch ${sampleId}_alternative_alleles_frequencies.png
    else
      predict_coinfection_illumina.py ${detected_variants_varscan_coinfection_txt} \
                                      ${sampleId} \
                                      /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS1.04_contamination.txt \
                                      /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS1.05_contamination.txt
    fi
    """
} 
