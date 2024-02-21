process lowCov {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(bam), path(bai)
    path ref_genome

    output:
    tuple val(sampleId), path('low_coverage.bed')
    tuple val(sampleId), path('lowcoverage_masked.fa')

    script:
    """
    pysam_quality_mask_final.py ${bam} 10 ${params.mask}
    cat quality_mask.bed | bedtools merge -d 2 | \
                           awk 'BEGIN {OFS = "\t"}; {if (\$3-\$2 > 3) print \$1,\$2,\$3}' >> low_coverage.bed
    bedtools maskfasta -fi ${ref_genome} \
                       -bed low_coverage.bed \
                       -fo lowcoverage_masked.fa
    """
}