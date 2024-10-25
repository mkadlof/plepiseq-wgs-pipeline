process masking {
    tag "masking:${sampleId}"

    input:
    tuple val(sampleId), path(bam), path(bai), path(primers), path(pairs), val(QC_status)

    output:
    tuple val(sampleId), path('ivar_trimmed_all.bam'), val(QC_status)

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch ivar_trimmed_all.bam
    else
      length=`echo "${params.length} - 40" | bc -l`
      ivar trim -i ${bam} \
               -b ${primers} \
               -m \${length} \
               -f ${pairs} \
               -q ${params.quality_initial} \
               -e \
               -p ivar_trimmed_all 
    fi
    """
}
