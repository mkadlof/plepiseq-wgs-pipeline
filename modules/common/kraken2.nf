process kraken2 {
    tag "kraken2:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}/kraken2", mode: 'copy', pattern: "report_kraken2_individualreads.txt"
    publishDir "${params.results_dir}/${sampleId}/kraken2", mode: 'copy', pattern: "report_kraken2.txt"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "summary_kraken.txt"
    containerOptions "--volume ${params.kraken2_db_absolute_path_on_host}:/home/external_databases/kraken2"
    maxForks 4

    input:
    tuple val(sampleId), path(reads), val(QC_STATUS)

    output:
    tuple val(sampleId), path('report_kraken2.txt'), path('report_kraken2_individualreads.txt'), path('Summary_kraken_genera.txt'), path('Summary_kraken_species.txt')

    script:
    """
    if [ ${QC_STATUS} == "nie" ]; then
        # upstream module failed we produce empty files do the pipeline can execute
        touch report_kraken2.txt
        touch report_kraken2_individualreads.txt
        touch Summary_kraken_genera.txt
        touch Summary_kraken_species.txt
    else
        kraken2 --db /home/external_databases/kraken2 \
                --report report_kraken2.txt \
                --threads ${params.threads} \
                --gzip-compressed \
                --minimum-base-quality ${params.quality_initial} \
                --use-names ${reads[0]} ${reads[1]} >> report_kraken2_individualreads.txt 2>&1

        # parse kraken extract two most abundant FAMILIES
        LEVEL="G" # G - genus, S - species
        SPEC1=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
        SPEC2=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
        ILE1==`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
        ILE2==`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

        echo -e "${sampleId}\t\${SPEC1}\${ILE1}%\t\${SPEC2}\${ILE2}%" >> Summary_kraken_genera.txt

        LEVEL="S"
        SPEC1=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
        SPEC2=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
        ILE1==`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
        ILE2==`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

        echo -e "${sampleId}\t\${SPEC1}\${ILE1}%\t\${SPEC2}\${ILE2}%" >> Summary_kraken_species.txt
    fi
    """
}
