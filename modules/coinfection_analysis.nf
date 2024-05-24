process coinfection_analysis {
    tag "Final step in coinfection analysis for sample:\t$sampleId"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy'

    input:
    tuple val(sampleId), path(detected_variants_varscan_contamination_txt)

    output:
    tuple val(sampleId), path("${sampleId}_alternative_alleles_frequencies.png"), path("${sampleId}_contamination_summary.txt")

    script:
    """
    # analiza konfekcji oparta jest o podobieństwo do 3 znanych sampli z EQA2023
    # nowa wersja skryptu pozawala podawac dowolna ilosc probek dla ktorych wiemy ze sa koinfekowane
    # wystarczy dodawac kolejen pliki
    predict_contamination_illumina.py ${detected_variants_varscan_contamination_txt} \
                                      ${sampleId} \
                                      /home/SARS-CoV2/contaminations/ESIB_EQA_2023.SARS2.09_contaminations.txt \
                                      /home/SARS-CoV2/contaminations/ESIB_EQA_2023.SARS2.17_contaminations.txt \
                                      /home/SARS-CoV2/contaminations/ESIB_EQA_2023.SARS2.32_contaminations.txt

    """
}
