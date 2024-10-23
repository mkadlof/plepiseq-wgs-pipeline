process kraken2_illumina {
    // This module combines old kraken2 with get_species_{illumina/nanopore} process from bacterial pipeline
    // It takes as an input reads and qc_status
    // and if qc status is not "nie" produces json (contaminations tab)
    // and updates qc_status (from "tak" to "nie") if provided criteria are not mean
    // criteria for now are extreamly liberal (at leat 5% of reads originate from a genus to which intendent species
    // belongs)
    tag "kraken2:${sampleId}"

    containerOptions "--volume ${params.kraken2_db_absolute_path_on_host}:/home/external_databases/kraken2"
    maxForks 4
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "Summary_kraken*"
    publishDir "${params.results_dir}/${sampleId}/json_output", mode: 'copy', pattern: "contaminations.json"
    input:
    tuple val(sampleId), path(reads), val(QC_STATUS)
    val(EXPECTED_GENUS)

    output:
    path('Summary_*'), emit: to_pubdir
    tuple val(sampleId), path('contaminations.json'), emit: json
    tuple val(sampleId), env(QC_status_contaminations), emit: qcstatus_only
    tuple val(sampleId), env(FINAL_GENUS), env(QC_status_contaminations), emit: species_and_qcstatus
    script:
    """
    if [ ${QC_STATUS} == "nie" ]; then
        # upstream module failed we produce empty files do the pipeline can execute
        QC_status_contaminations="nie"
        FINAL_GENUS="unknown"
        touch  report_kraken2.txt
        ERR_MSG="This module recieved wrong QC status and did not produce any valid output"
        json_output_contaminations.py -k report_kraken2.txt -g skip -x skip -y skip -s nie -m "\${ERR_MSG}" -o contaminations.json

    else
        kraken2 --db /home/external_databases/kraken2 \
                --report report_kraken2.txt \
                --threads ${params.threads} \
                --gzip-compressed \
                --minimum-base-quality ${params.quality_initial} \
                --use-names ${reads[0]} ${reads[1]} >> report_kraken2_individualreads.txt 2>&1

        # parse kraken extract two most abundant Genus and Ssecies for json output
        # This
        LEVEL="G" # G - genus, S - species
        GENUS1=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
        GENUS2=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
        GENUS1_ILE=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
        GENUS2_ILE=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

        echo -e "${sampleId}\t\${GENUS1}\${GENUS1_ILE}%\t\${GENUS2}\${GENUS2_ILE}%" >> Summary_kraken_genera.txt

        LEVEL="S"
        SPEC1=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
        SPEC2=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
        SPEC1_ILE=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
        SPEC2_ILE==`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

        echo -e "${sampleId}\t\${SPEC1}\${SPEC1_ILE}%\t\${SPEC2}\${SPEC2_ILE}%" >> Summary_kraken_species.txt
        # Extract number of reads associated with expected genus to check if the number of reads is greater than params
        # and return failed QC if sample does not meet this criterion
        GENUS_EXPECTED_ILE=`cat report_kraken2.txt | grep -w G |  grep "${EXPECTED_GENUS}" | tr -s " " | cut -f1 | tr -d " " | awk '{print int(\$1)}'`

        # This section introducec criteria to switch from "tak" to "nie" in this module
        # all downstream modules will not execute and produce dummy values
        if [ \${GENUS_EXPECTED_ILE} -lt ${params.expected_genus_value} ]; then
            QC_status_contaminations="nie"
            FINAL_GENUS="${EXPECTED_GENUS}"
            ERR_MSG="Number of reads associated with the expected genus is below threshold set to: ${params.expected_genus_value} %"
            json_output_contaminations.py -k report_kraken2.txt -g skip -x skip -y skip -s blad -m "\${ERR_MSG}" -o contaminations.json

        else
            FINAL_GENUS="unknown"
            QC_status_contaminations="tak"
            json_output_contaminations.py -k report_kraken2.txt -g skip -x skip -y skip -s tak -o contaminations.json

        fi


    fi
    """
}

process kraken2_nanopore {
    // This process only differ from illumina in kraken2 execution
    tag "kraken2:${sampleId}"

    containerOptions "--volume ${params.kraken2_db_absolute_path_on_host}:/home/external_databases/kraken2"
    maxForks 4
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "Summary_kraken*"
    publishDir "${params.results_dir}/${sampleId}/json_output", mode: 'copy', pattern: "contaminations.json"
    input:
    tuple val(sampleId), path(reads), val(QC_STATUS)
    val(EXPECTED_GENUS)

    output:
    path('Summary_*'), emit: to_pubdir
    tuple val(sampleId), path('contaminations.json'), emit: json
    tuple val(sampleId), env(QC_status_contaminations), emit: qcstatus_only
    tuple val(sampleId), env(FINAL_GENUS), env(QC_status_contaminations), emit: species_and_qcstatus
    script:
    """
    if [ ${QC_STATUS} == "nie" ]; then
        # upstream module failed we produce empty files do the pipeline can execute
        QC_status_contaminations="nie"
        FINAL_GENUS="unknown"
        touch  report_kraken2.txt
        ERR_MSG="This module recieved wrong QC status and did not produce any valid output"
        json_output_contaminations.py -k report_kraken2.txt -g skip -x skip -y skip -s nie -m "\${ERR_MSG}" -o contaminations.json

    else
        kraken2 --db /home/external_databases/kraken2 \
                --report report_kraken2.txt \
                --threads ${params.threads} \
                --gzip-compressed \
                --minimum-base-quality ${params.quality_initial} \
                --use-names ${reads} >> report_kraken2_individualreads.txt 2>&1

        # parse kraken extract two most abundant Genus and Ssecies for json output
        # This
        LEVEL="G" # G - genus, S - species
        GENUS1=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
        GENUS2=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
        GENUS1_ILE=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
        GENUS2_ILE=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

        echo -e "${sampleId}\t\${GENUS1}\${GENUS1_ILE}%\t\${GENUS2}\${GENUS2_ILE}%" >> Summary_kraken_genera.txt

        LEVEL="S"
        SPEC1=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
        SPEC2=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
        SPEC1_ILE=`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
        SPEC2_ILE==`cat report_kraken2.txt | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

        echo -e "${sampleId}\t\${SPEC1}\${SPEC1_ILE}%\t\${SPEC2}\${SPEC2_ILE}%" >> Summary_kraken_species.txt
        # Extract number of reads associated with expected genus to check if the number of reads is greater than params
        # and return failed QC if sample does not meet this criterion
        GENUS_EXPECTED_ILE=`cat report_kraken2.txt | grep -w G |  grep "${EXPECTED_GENUS}" | tr -s " " | cut -f1 | tr -d " " | awk '{print int(\$1)}'`

        # This section introducec criteria to switch from "tak" to "nie" in this module
        # all downstream modules will not execute and produce dummy values
        if [ \${GENUS_EXPECTED_ILE} -lt ${params.expected_genus_value} ]; then
            QC_status_contaminations="nie"
            FINAL_GENUS="${EXPECTED_GENUS}"
            ERR_MSG="Number of reads associated with the expected genus is below threshold set to: ${params.expected_genus_value} %"
            json_output_contaminations.py -k report_kraken2.txt -g skip -x skip -y skip -s blad -m "\${ERR_MSG}" -o contaminations.json

        else
            FINAL_GENUS="unknown"
            QC_status_contaminations="tak"
            json_output_contaminations.py -k report_kraken2.txt -g skip -x skip -y skip -s tak -o contaminations.json

        fi


    fi
    """
}
