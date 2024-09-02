process lowCov {
    tag "lowCov:${sampleId}"

    input:
    tuple val(sampleId), path(bam), path(bai), path(ref_genome)

    output:
    tuple val(sampleId), path('low_coverage.bed')
    tuple val(sampleId), path('lowcoverage_masked.fa')

    script:
    """
    position_quality_for_coverage=10
    predict_lowcoverage_pysam.py ${bam} \${position_quality_for_coverage} ${params.mask} ${ref_genome}
    for K in `ls *mask.bed`; do cat \${K} | bedtools merge -d 2 | awk 'BEGIN {OFS = "\t"}; {if (\$3-\$2 >= 2) print \$1,\$2,\$3}' >> low_coverage.bed; done
    bedtools maskfasta -fi ${ref_genome} \
                       -bed low_coverage.bed \
                       -fo lowcoverage_masked.fa
    """
}
