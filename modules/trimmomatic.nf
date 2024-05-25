process trimmomatic {
    tag "Initial filtering nad clipping of reads for sample:\t$sampleId"
    //publishDir "${params.results_dir}/${sampleId}/reads_clipping", mode: 'symlink'

    input:
    tuple val(sampleId), path(reads)
    val(adapters)

    output:
    tuple val(sampleId), path('*_paired.fastq.gz')
    tuple val(sampleId), path('*_unpaired.fastq.gz')

    script:
    """
    java -jar /opt/trimmomatic/trimmomatic.jar PE ${reads[0]} ${reads[1]} \
                                            forward_paired.fastq.gz \
                                            forward_unpaired.fastq.gz \
                                            reverse_paired.fastq.gz \
                                            reverse_unpaired.fastq.gz \
                                            ILLUMINACLIP:${adapters}:2:30:10:8:True \
                                            LEADING:${params.quality_initial} \
                                            TRAILING:${params.quality_initial} \
                                            SLIDINGWINDOW:4:4 \
                                            MINLEN:${params.length}
    """
}
