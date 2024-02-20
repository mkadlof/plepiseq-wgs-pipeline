process trimmomatic {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(reads)
    path adapters

    output:
    tuple val(sampleId), path('forward_paired.fq.gz'), path('forward_unpaired.fq.gz'), path('reverse_paired.fq.gz'), path('reverse_unpaired.fq.gz')

    script:
    """
    java -jar /opt/trimmomatic/trimmomatic.jar PE ${reads[0]} ${reads[1]} \
                                            forward_paired.fq.gz \
                                            forward_unpaired.fq.gz \
                                            reverse_paired.fq.gz \
                                            reverse_unpaired.fq.gz \
                                            ILLUMINACLIP:${adapters}:2:30:10:8:True \
                                            LEADING:${params.quality_initial} \
                                            TRAILING:${params.quality_initial} \
                                            SLIDINGWINDOW:4:4 \
                                            MINLEN:${params.length}
    """
}