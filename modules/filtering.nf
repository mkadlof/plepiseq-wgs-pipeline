process filtering {
    publishDir "results/${sampleId}", mode: 'symlink', pattern: 'Statystyki_one_amplicon.txt'

    input:
    tuple val(sampleId), path(bam), path(bai)
    path primers

    output:
    tuple val(sampleId), path('reads_inneramplicon.bam'), path('reads_inneramplicon.bam.bai'),
                         path('clean_sort_dedup_trimmed_sort.bam'), path('clean_sort_dedup_trimmed_sort.bam.bai'),
                         path('Statystyki_one_amplicon.txt')

    script:
    """
    simple_filter_amplicon_10_illumina.py ${bam} ${primers}
    samtools sort -@ ${params.threads} -o reads_inneramplicon.bam reads_inneramplicon.bam
    samtools index reads_inneramplicon.bam

    TO_MERGE=`ls -l reads_two_amplicons_for_*sorted_ivar.bam 2> /dev/null | tr -s " " |  cut -d " " -f9 | tr "\n" " "`

    if [ -n "\${TO_MERGE}" ]; then
        samtools merge -o clean_sort_dedup_trimmed_sort.bam \${TO_MERGE}
    else
        samtools view -b -o clean_sort_dedup_trimmed_sort.bam -L /dev/null reads_inneramplicon.bam
    fi

    samtools sort -@ ${params.threads} -o clean_sort_dedup_trimmed_sort.bam clean_sort_dedup_trimmed_sort.bam
    samtools index clean_sort_dedup_trimmed_sort.bam
    """
}