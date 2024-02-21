process filtering {
    publishDir "results/${sampleId}", mode: 'symlink', pattern: 'Statystyki_one_amplicon.txt'

    input:
    tuple val(sampleId), path(bam), path(bai)
    path primers

    output:
    tuple val(sampleId), path('reads_inneramplicon.bam'), path('reads_inneramplicon.bam.bai')
    tuple val(sampleId), path('tmp.bam') //, path('tmp.bam.bai')
    tuple val(sampleId), path('Statystyki_one_amplicon.txt')

    script:
    """
    simple_filter_amplicon_10_illumina.py ${bam} ${primers}
    samtools sort -@ ${params.threads} -o reads_inneramplicon.bam reads_inneramplicon.bam
    samtools index reads_inneramplicon.bam

    TO_MERGE=`ls -l reads_two_amplicons_for_*sorted_ivar.bam 2> /dev/null | tr -s " " |  cut -d " " -f9 | tr "\n" " "`

    if [ -n "\${TO_MERGE}" ]; then
        samtools merge -o tmp.bam \${TO_MERGE}
    else
        samtools view -b -o tmp.bam -L /dev/null reads_inneramplicon.bam
    fi

    samtools sort -@ ${params.threads} -o tmp.bam tmp.bam
    # samtools index tmp.bam
    """
}