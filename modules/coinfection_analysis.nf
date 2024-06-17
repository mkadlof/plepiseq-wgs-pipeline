process coinfection_analysis {
    tag "coinfection_analysis:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy'

    input:
    tuple val(sampleId), path(detected_variants_varscan_contamination_txt)

    output:
    tuple val(sampleId), path("${sampleId}_alternative_alleles_frequencies.png"), path("${sampleId}_contamination_summary.txt")

    script:
    """
    # The coinfection analysis is based on similarity to 3 known samples from EQA2023.
    # The new version of the script allows specifying any number of samples that are known to be
    # coinfected. Simply add additional files.

    predict_contamination_illumina.py ${detected_variants_varscan_contamination_txt} \
                                      ${sampleId} \
                                      /home/data/sarscov2/contaminations/ESIB_EQA_2023.SARS2.09_contaminations.txt \
                                      /home/data/sarscov2/contaminations/ESIB_EQA_2023.SARS2.17_contaminations.txt \
                                      /home/data/sarscov2/contaminations/ESIB_EQA_2023.SARS2.32_contaminations.txt

    """
}
