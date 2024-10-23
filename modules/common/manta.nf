process introduce_SV_with_manta {
    tag "manta:$sampleId"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy'
    container = params.manta_image

    input:
    tuple val(sampleId), path(bam_file), path(bai_file),  path(ref_genome_with_index), val(QC_status_picard), path(consensus_masked_fasta), val(QC_status_consensus)

    output:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path(ref_genome_with_index), env(QC_status_exit), emit: fasta_refgenome_and_qc
    tuple val(sampleId), path('dummy_json.json'), emit: json

    script:
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }

    """
    if [[ ${QC_status_picard} == "nie"  || ${QC_status_consensus} == "nie" ]]; then
      touch output_consensus_masked_SV.fa
      touch dummy_json.json
      QC_status_exit="nie"
    else
      touch dummy_json.json
      QC_status_exit="tak"
      samtools faidx  ${ref_genome_with_index[final_index]}
      python /opt/docker/manta/bin/configManta.py --bam ${bam_file} --reference ${ref_genome_with_index[final_index]} --runDir Manta_results
      python Manta_results/runWorkflow.py -j ${params.threads} --quiet

      if [ -e Manta_results/results/variants/diploidSV.vcf.gz ]; then
          # Wywalamy skomplikowane SV jak translokacje itd typun BND
          bcftools view -O z -o manta_results.vcf.gz -i '(FILTER="PASS" | FILTER="MaxDepth" | FILTER="NoPairSupport") && SVTYPE != "BND"' Manta_results/results/variants/diploidSV.vcf.gz
          tabix manta_results.vcf.gz

          ILE_SV=`zcat manta_results.vcf.gz  | grep SVTYPE | grep -v INFO | wc -l`

          if [ \${ILE_SV} -gt 0 ]; then
              
              cat ${ref_genome_with_index[final_index]} | bcftools consensus -s - manta_results.vcf.gz  > output_manta.fa
              /home/bin/sarscov2/insert_SV_python2.py ${consensus_masked_fasta} output_manta.fa output_consensus_masked_SV.fa
          else
              HEADER=`head -1 ${consensus_masked_fasta}`
              NEW_HEADER=`echo -e "\${HEADER}_SV"`
              cat ${consensus_masked_fasta} | sed s"/\${HEADER}/\${NEW_HEADER}/"g > output_consensus_masked_SV.fa

          fi
      else
        QC_status_exit="nie"
        touch output_consensus_masked_SV.fa
      fi
    fi
    """
}
