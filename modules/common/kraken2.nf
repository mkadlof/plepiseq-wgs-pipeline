process kraken2 {
    tag "kraken2:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}/kraken2", mode: 'copy', pattern: "report_kraken2_individualreads.txt"
    publishDir "${params.results_dir}/${sampleId}/kraken2", mode: 'copy', pattern: "report_kraken2.txt"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "summary_kraken.txt"
    containerOptions "--volume ${params.kraken2_db_absolute_path_on_host}:/home/external_databases/kraken2"
    maxForks 5

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path('report_kraken2.txt'), path('report_kraken2_individualreads.txt')

    script:
    """
    kraken2 --db /home/external_databases/kraken2 \
            --report report_kraken2.txt \
            --threads ${params.threads} \
            --gzip-compressed \
            --minimum-base-quality ${params.quality_initial} \
            --use-names ${reads[0]} ${reads[1]} >> report_kraken2_individualreads.txt 2>&1

    # parse kraken extract two most abundant species
    SPEC1=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w S | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
    SPEC2=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w S | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`

    ILE1==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w S | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
    ILE2==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w S | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

    echo -e "${sampleId}\t\${SPEC1}\${ILE1}%\t\${SPEC2}\${ILE2}%" >> summary_kraken.txt
    """
}
