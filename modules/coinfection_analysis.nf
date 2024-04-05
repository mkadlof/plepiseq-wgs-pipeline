process coinfection_analysis {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(detected_variants_varscan_contamination_txt)

    output:
    tuple val(sampleId), path("${sampleId}_alternative_alleles_frequencies.png"), path("${sampleId}_contamination_summary.txt")

    script:
    """
    # analiza konfekcji oparta jest o podobieństwo do 3 znanych sampli z EQA2023
    predict_contamination_illumina.py ${detected_variants_varscan_contamination_txt} \
                                      ${sampleId} \
                                      /home/SARS-CoV2/contaminations/ESIB_EQA_2023.SARS2.09_contaminations.txt \
                                      /home/SARS-CoV2/contaminations/ESIB_EQA_2023.SARS2.17_contaminations.txt \
                                      /home/SARS-CoV2/contaminations/ESIB_EQA_2023.SARS2.32_contaminations.txt

    """
}