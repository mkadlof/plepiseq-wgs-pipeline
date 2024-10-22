process picard_downsample {
    tag "picard:${sampleId}"

    input:
    tuple val(sampleId), path(bam), path(bai), val(QC_status)

    output:
    tuple val(sampleId), path('downsample.bam'), path('downsample.bam.bai'), env(QC_exit), emit: to_pubdir

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch downsample.bam
      touch downsample.bam.bai
      QC_exit="nie"
    else

      ILE=`samtools view ${bam} | wc -l `
      ILE=`echo "\${ILE} + 2" | bc -l`
      NORM=`echo "${params.max_number_for_SV}/\${ILE}" | bc -l `
      if [ `awk -v norm="\${NORM}" 'BEGIN {if(norm>0.99) {print 1} else {print 0}}'` -eq 1 ]; then 
        NORM=0.99
      fi

      suffix=`echo "${params.max_number_for_SV}/1000" | bc -l | awk '{print int(\$0)}'`

      java -jar /opt/picard/picard.jar PositionBasedDownsampleSam \
                                      --INPUT ${bam} \
                                      --OUTPUT downsample.bam -F \${NORM}
      samtools index downsample.bam
      NO_READS=`samtools view downsample.bam | wc -l`
      if [ "\${NO_READS}" -lt 1000 ]; then
        QC_exit="nie"
      else
        QC_exit="tak"
      fi # koniec if na za malo odczytow w downsampled bam
    
    fi # koniec if na initial QC_status
    """
}
