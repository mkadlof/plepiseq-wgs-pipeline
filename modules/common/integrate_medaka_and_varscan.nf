process medaka_varscan_integration {
    // This module integrates medaka vcf and SPECIFIC parts of varscan prediction
    // i.e. SNPs with over 80% usage (according to varscan) that are not predicted as valid SNPs according to medaka
    // due to complex underlying genotype of a region
    tag "programs integration:${sampleId}"
    container  = params.main_image

    input:
    tuple val(sampleId), path('medaka_annotated_filtered.vcf.gz'), path('medaka_annotated_filtered.vcf.gz.tbi'), path('genome.fasta'), val(QC_status),  path('detected_variants_varscan.txt')

    output:
    tuple val(sampleId), path('medaka_and_varscan_final.vcf.gz'), path('medaka_and_varscan_final.vcf.gz.tbi'), val(QC_status), emit: vcf
    tuple val(sampleId),  path('genome.fasta'),  val(QC_status), emit: reference_genome

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch medaka_and_varscan_final.vcf.gz
      touch medaka_and_varscan_final.vcf.gz.tbi
    else
        MINIMUM_USAGE=80
        awk -v usage="\${MINIMUM_USAGE}" '{if (int(\$7) > usage) print \$0; else if (\$1 ~ /Chrom/) print \$0}' detected_variants_varscan.txt >> detected_variants_varscan_part1_filtered.txt
        
        merge_varscan_with_medaka_final_INFL_first_round.py medaka_annotated_filtered.vcf.gz detected_variants_varscan_part1_filtered.txt  medaka_and_varscan.vcf

        bcftools sort medaka_and_varscan.vcf |  bcftools norm -c w -d all -f genome.fasta | bcftools norm -c w -m -indels -f genome.fasta | bcftools filter -O z -o medaka_and_varscan_final.vcf.gz -i "QUAL >= 0 && INFO/DP >= 1"

        tabix medaka_and_varscan_final.vcf.gz
    


    fi
    """
}
