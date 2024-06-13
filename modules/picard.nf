process picard {
    tag "picard:${sampleId}"

    input:
    tuple val(sampleId), path(bam), path(bai)
    val(primers)
    val(pairs)

    output:
    tuple val(sampleId), path('downsample.bam'), path('downsample.bam.bai')

    script:
    """
    ILE=`samtools view ${bam} | wc -l `
    # dodajemy 2 żeby nie miec pustych przelotów
    ILE=`echo "\${ILE} + 2" | bc -l`
    NORM=`echo "${params.max_number_for_SV}/\${ILE}" | bc -l `
    if (( \$(echo "\${NORM} > 0.99" | bc -l ) )); then
            NORM=0.99
    fi

    suffix=`echo "${params.max_number_for_SV}/1000" | bc -l | awk '{print int(\$0)}'`

    java -jar /opt/picard/picard.jar PositionBasedDownsampleSam \
                                    --INPUT ${bam} \
                                    --OUTPUT downsample.bam -F \${NORM}
    samtools index downsample.bam
    """
}
