process viterbi {
    tag "Pre-lofreq steps for sample:\t${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "forvariants.bam*"

    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(sampleId), path('forvariants.bam'), path('forvariants.bam.bai')

    script:
    """
    #lofreq viterbi --ref ${reference_fasta} \
    #               --out clean_sort_dedup_trimmed_sort_viterbi.bam \
    #               ${bam}

    #lofreq indelqual --ref ${reference_fasta} \
    #                 --out forvariants.bam \
    #                 --dindel clean_sort_dedup_trimmed_sort_viterbi.bam

    # Zmiana poniewaz w sample 05 w poblizu regionu 22000 gdzie mamy delecje blisko 5' odczytu  + mutacje tuz obok
    # Program tak "przesuwal" alignment ze "zasypal" czesciowo delecje prowadzac do braku identyfikacji delecji 
    # przez lofreq i varscan
   
    lofreq indelqual --ref ${reference_fasta} \
                     --out forvariants.bam \
                     --dindel ${bam}

    samtools sort -@ ${params.threads} \
                  -o forvariants.bam \
                  forvariants.bam
    samtools index forvariants.bam
    """
}
