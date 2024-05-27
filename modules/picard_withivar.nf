process picard {
    tag "Preparing bam for SV caller for sample:\t$sampleId"
    // publishDir "${params.results_dir}/${sampleId}/bam_for_SV", mode: 'symlink', pattern: "downsample.bam*"
    input:
    tuple val(sampleId), path(bam), path(bai)
    val(primers)
    val(pairs)
    output:
    tuple val(sampleId), path('downsample_ivar_sorted.bam'), path('downsample_ivar_sorted.bam.bai')

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


    # Dodajemy rowniez usuuwanie primerow w tej sekwencji
    # ALe musimy obejsc sytuacje w ktorej nie ma mapowan
    

    if [ \${ILE} -lt 100 ] ;then
        cp downsample.bam downsample_ivar_sorted.bam
        cp downsample.bam.bai downsample_ivar_sorted.bam.bai
    else
    ivar trim -i downsample.bam \
             -b ${primers} \
             -q ${params.quality_initial} \
             -e \
             -p downsample_ivar
    
    samtools sort -@ ${params.threads} -o downsample_ivar_sorted.bam downsample_ivar.bam
    samtools index downsample_ivar_sorted.bam
    fi
    # MAnta nie rozumie chyba soft-clipping 
    # patrz przykladi 06 z SARS2 z EQA2024
    # z tego powodu wracam do prostszej wersji modulu
    """
}
