# Statistics: Modules wgsMetrics and simpleStats


## Module wgsMetrics
```Bash
    java -jar /opt/picard/picard.jar CollectWgsMetrics --REFERENCE_SEQUENCE ${reference_fasta} \
                                                   --MINIMUM_BASE_QUALITY ${params.quality_initial} \
                                                   --MINIMUM_MAPPING_QUALITY 30 \
                                                   --INPUT ${bam} \
                                                   --OUTPUT picard_statistics.txt
```

## Module simpleStats

```Bash
   calculate_N.py ${consensus_masked_fa} consensus_masked_N_summary.txt
    
    # wybrane pozycje z picard-a
    OUT=`cat ${picard_statistics_txt} | head -8 | tail -1`
    HEADER=`cat ${picard_statistics_txt} | head -7 | tail -1`
    
    echo -e "\${HEADER}\t\${OUT}" > picard_summary.txt
    
    # użycie primerów
    touch log.txt
    MES1=`primer_usage_sum.py ${params.primers} log.txt 40 | head -1`
    MES2=`primer_usage_sum.py ${params.primers} log.txt 40 | head -2 | tail -1`
    
    echo -e "\${MES1}" > primers_poor_stretch.txt
    echo -e "\${MES2}" > primers_poor.txt
```

Here we parse partial results/logs from programs:
1. Counting the number of N's in the genome.
2. Summary of the percentage of the genome covered at 5x, 10x, 20x, 30x.
3. Summary of which primers had low usage.
