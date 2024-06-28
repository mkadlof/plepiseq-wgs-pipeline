process masking {
    tag "masking:${sampleId1}"

    input:
    tuple val(sampleId1), path(bam), path(bai)
    tuple val(sampleId2), path(primers_and_pairs)

    output:
    tuple val(sampleId1), path('ivar_trimmed_all.bam')

    script:
    """
    length=`echo "${params.length} - 40" | bc -l`
    ivar trim -i ${bam} \
             -b ${primers_and_pairs[0]} \
             -m \${length} \
             -f ${primers_and_pairs[1]} \
             -q ${params.quality_initial} \
             -e \
             -p ivar_trimmed_all
    """
}
