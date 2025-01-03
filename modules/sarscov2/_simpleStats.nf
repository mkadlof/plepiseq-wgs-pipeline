process simpleStats {
    tag "simpleStats:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "primers_poor*"

    input:
    tuple val(sampleId), path(consensus_masked_fa), path(picard_statistics_txt), path(primers)

    output:
    tuple val(sampleId), path('consensus_masked_N_summary.txt'), path('picard_summary.txt'), path('primers_poor_stretch.txt'), path('primers_poor.txt')

    script:
    """
    calculate_N.py ${consensus_masked_fa} consensus_masked_N_summary.txt
    
    # wybrane pozycje z picard-a
    OUT=`cat ${picard_statistics_txt} | head -8 | tail -1`
    HEADER=`cat ${picard_statistics_txt} | head -7 | tail -1`
    
    echo -e "\${HEADER}\t\${OUT}" > picard_summary.txt
    
    # użycie primerów
    touch log.txt
    echo ${primers}
    MES1=`primer_usage_sum.py ${primers} log.txt 40 | head -1`
    MES2=`primer_usage_sum.py ${primers} log.txt 40 | head -2 | tail -1`
    
    echo -e "\${MES1}" > primers_poor_stretch.txt
    echo -e "\${MES2}" > primers_poor.txt
    """
}
