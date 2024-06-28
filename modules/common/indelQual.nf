process indelQual {
    tag "indelQual:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "forvariants.bam*"

    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple val(sampleId2), path(ref_genome)

    output:
    tuple val(sampleId), path('forvariants.bam'), path('forvariants.bam.bai')

    script:
    """
    lofreq indelqual --ref ${ref_genome} \
                     --out forvariants.bam \
                     --dindel ${bam}
    samtools sort -@ ${params.threads} \
                  -o forvariants.bam \
                  forvariants.bam
    samtools index forvariants.bam
    """
}
