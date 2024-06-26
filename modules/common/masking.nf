process masking {
    tag "masking:${sampleId1}"

    input:
    tuple val(sampleId1), path(bam), path(bai)
    tuple val(sampleId2), path(primers)
    tuple val(sampleId3), path(pairs)

    output:
    tuple val(sampleId1), path('ivar_trimmed_all.bam')

    script:
    """
    length=`echo "${params.length} - 40" | bc -l`
    ivar trim -i ${bam} \
             -b ${primers} \
             -m \${length} \
             -f ${pairs} \
             -q ${params.quality_initial} \
             -e \
             -p ivar_trimmed_all
    """
}
