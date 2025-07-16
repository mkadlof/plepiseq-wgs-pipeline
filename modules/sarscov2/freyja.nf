process freyja_sars {
    tag "freyja:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    // publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "coinfections.tsv"
    containerOptions "--volume ${params.external_databases_path}:/home/external_databases/"

    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), path(ref_genome_with_index), val(QC_status)

    output:
    // tuple val(sampleId), path('coinfections.tsv'), emit: to_pubdir
    tuple val(sampleId), path('coinfections_freyja.json'), emit: json

    script:
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }
    """
    if [ ${QC_status} == "nie" ]; then
      # echo "No reads in the bam file"
      touch coinfections.tsv
      freyja_status="nie"
      if [ "${params.lan}" == "pl" ]; then
        ERR_MSG="Ten moduł został uruchomiony na próbce, która nie przeszła kontroli jakości."
      else
        ERR_MSG="This sample failed a QC analysis during an earlier phase of the analysis."
      fi
      
      echo -e "{\\"status\\":\\"\${freyja_status}\\",
                \\"error_message\\":\\"\${ERR_MSG}\\"}" >> coinfections_freyja.json
    else
      mkdir variants_files depth_files demix_files
      MEAN_DEPTHS=`samtools depth -aa mapped_reads.bam | awk '{sum+=\$3} END { print int(sum/NR)}'`
      # We call variants with freyja only if mean depth is above 20
      if [ \${MEAN_DEPTHS} -gt 20 ]; then


        freyja variants mapped_reads.bam --minq ${params.freyja_minq} --variants variants_files/test.variants.tsv --depths depth_files/test.depth --ref ${ref_genome_with_index[final_index]}
        freyja demix variants_files/test.variants.tsv depth_files/test.depth --output demix_files/test.output --confirmedonly --barcodes  /home/external_databases/freyja/sarscov2/usher_barcodes.csv
        freyja aggregate demix_files/ --output coinfections.tsv
        freyja_status="tak"
        freyja_lineage_1_name=`cat coinfections.tsv  | cut -f3 | tail -1 | cut -d " " -f1`
        freyja_lineage_2_name=`cat coinfections.tsv  | cut -f3 | tail -1 | cut -d " " -f2`
        freyja_lineage_1_abundance=`cat coinfections.tsv  | cut -f4 | tail -1 | cut -d " " -f1 | awk '{printf "%.2f", \$0}'`
        freyja_lineage_2_abundance=`cat coinfections.tsv  | cut -f4 | tail -1 | cut -d " " -f2 | awk '{printf "%.2f", \$0}'`

        # In case freyja found a single linage
        if [ -z \${freyja_lineage_2_name} ]; then
          freyja_lineage_2_name="unk"
          freyja_lineage_2_abundance=0
        fi

        # In case both linages have the same name
        if [ \${freyja_lineage_1_name} == \${freyja_lineage_2_name} ]; then
          freyja_lineage_2_name="unk"
          freyja_lineage_2_abundance=0
        fi

        echo -e "{\\"status\\":\\"\${freyja_status}\\",
                  \\"freyja_lineage1_name\\":\\"\${freyja_lineage_1_name}\\",
                  \\"freyja_lineage2_name\\":\\"\${freyja_lineage_2_name}\\",
                  \\"freyja_lineage1_abundance\\":\${freyja_lineage_1_abundance},
                  \\"freyja_lineage2_abundance\\":\${freyja_lineage_2_abundance}}" >> coinfections_freyja.json
      else
        if [ "${params.lan}" == "pl" ]; then
          ERR_MSG="Średnie pokrycie dla tej próbki: \${MEAN_DEPTHS}, jest poniżej zalecanego progu przez Freya"
        else
           ERR_MSG="Mean depth for this sample: \${MEAN_DEPTHS}, which is below Freyja recommended threshold"
         fi
           touch coinfections.tsv
           freyja_status="nie"
           echo -e "{\\"status\\":\\"\${freyja_status}\\",
                     \\"error_message\\":\\"\${ERR_MSG}\\"}" >> coinfections_freyja.json
        
      fi
    fi

    """
}
