process filtering {
    tag "filtering:${sampleId}"

    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple val(sampleId2), path(primers)

    output:
    tuple val(sampleId), path('first_pass_sorted.bam'), path('first_pass_sorted.bam.bai') 
    tuple val(sampleId), path('two_amplicons_sorted.bam')
    tuple val(sampleId), path('Statystyki.txt')

    script:
    """
    bed_offset=0
    length=`echo "${params.length} - 20" | bc -l`
    equal_depth=`echo "${params.max_depth} / 2" | bc -l | awk '{print int(\$0)}'` 
    simple_filter_illumina_one_segment.py ${bam} ${primers} \${bed_offset} \${length} ${params.min_mapq} ${params.window_size} \${equal_depth}
    samtools index first_pass_sorted.bam
    
    if [ -e two_amplicons_sorted.bam ]; then
        cp two_amplicons_sorted.bam tmp.bam
    else
        samtools view -b -o two_amplicons_sorted.bam -L /dev/null first_pass_sorted.bam
    fi
    """
}
