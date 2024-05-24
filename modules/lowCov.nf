process lowCov {
    tag "Predicting low coverage regions for sample:\t$sampleId"
    //publishDir "${params.results_dir}/${sampleId}/lowcoverage", mode: 'symlink'

    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(sampleId), path('low_coverage.bed')
    tuple val(sampleId), path('lowcoverage_masked.fa')

    script:
    """
    position_quality_for_coverage=10
    predict_lowcoverage_pysam.py ${bam} \${position_quality_for_coverage} ${params.mask} ${reference_fasta}
    cat quality_mask.bed | bedtools merge -d 2 | \
                           awk 'BEGIN {OFS = "\t"}; {if (\$3-\$2 > 3) print \$1,\$2,\$3}' >> low_coverage.bed
    bedtools maskfasta -fi ${reference_fasta} \
                       -bed low_coverage.bed \
                       -fo lowcoverage_masked.fa
    """
}
