process manta {
    tag "manta:$sampleId"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy'
    container 'nf_illumina_sars-3.0-manta:latest'

    input:
    tuple val(sampleId), path(bam_file), path(bai_file), path(consensus_masked_fasta)

    output:
    tuple val(sampleId), path('output_consensus_masked_SV.fa')

    script:
    """
    GENOME_ID="MN908947.3"
    GENOME_FASTA="/SARS-CoV2/genome/sarscov2.fasta"
    ILE_ODCZYTOW=`samtools view ${bam_file} | wc -l`
    if [  \${ILE_ODCZYTOW} -lt 1000 ]; then
        # pusty bam, nie puszczamy manty, po prostu tworzymy kopie plikow z poprawionymi nazwami
        HEADER=`head -1 ${consensus_masked_fasta}`
        NEW_HEADER=`echo -e "\${HEADER}_SV"`
        cat ${consensus_masked_fasta} | sed s"/\${HEADER}/\${NEW_HEADER}/"g > output_consensus_masked_SV.fa

        #chmod a+wr output_consensus_masked_SV.fa

    else
        python /opt/docker/manta/bin/configManta.py --bam ${bam_file} --reference \${GENOME_FASTA} --runDir Manta_results
        python Manta_results/runWorkflow.py -j ${params.threads} --quiet

        if [ -e Manta_results/results/variants/diploidSV.vcf.gz ]; then
            # Wywalamy skomplikowane SV jak translokacje itd typun BND
            bcftools view -O z -o manta_results.vcf.gz -i '(FILTER="PASS" | FILTER="MaxDepth" | FILTER="NoPairSupport") && SVTYPE != "BND"' Manta_results/results/variants/diploidSV.vcf.gz
                    tabix manta_results.vcf.gz

            ILE_SV=`zcat manta_results.vcf.gz  | grep \${GENOME_ID} | grep -v cont | wc -l`

            if [ \${ILE_SV} -gt 0 ]; then
                cat \${GENOME_FASTA} | bcftools consensus -s - manta_results.vcf.gz  > output_manta.fa
                /home/bin/sarscov2/insert_SV_python2.py ${consensus_masked_fasta} output_manta.fa output_consensus_masked_SV.fa
            else
                HEADER=`head -1 ${consensus_masked_fasta}`
                NEW_HEADER=`echo -e "\${HEADER}_SV"`
                cat ${consensus_masked_fasta} | sed s"/\${HEADER}/\${NEW_HEADER}/"g > output_consensus_masked_SV.fa

            fi
        else
            echo "Error brak pliku Manta_results/results/variants/diploidSV.vcf.gz"
            exit 1
        fi
    fi
    """
}
