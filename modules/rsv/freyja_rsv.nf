process freyja_rsv {
    // We need to repeat mapping for RSV as Freyja needs a different genome than we used in main pipeline
    // We also "attach" this module not AFTER bwa/minmap2 but NEXT to it
    // We do not mask primers here
 
    tag "freyja:${sampleId}"
    container  = params.main_image
    cpus params.threads
    memory "40 GB"
    containerOptions "--volume ${params.external_databases_path}:/home/external_databases/"

    input:
    tuple val(sampleId), path(reads), val(TYPE), val(QC_status)

    output:
    tuple val(sampleId), path('coinfections_freyja.json'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
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

      cp /home/external_databases/freyja/RSV_${TYPE}/reference.fasta .
      bwa index reference.fasta

      # We map  reads all the reads to Freyja-required genome 
      if [ ${params.machine} == 'Illumina' ]; then  
        bwa mem -t ${task.cpus} -T 30 reference.fasta ${reads[0]} ${reads[1]} | \
        samtools view -@ ${task.cpus} -Sb -f 3 -F 2048 - | \
        samtools sort -@ ${task.cpus} -o mapped_reads.bam -
        samtools index mapped_reads.bam
      elif [ ${params.machine} == 'Nanopore' ]; then
        minimap2 -a -x map-ont -t ${task.cpus} -o tmp.sam reference.fasta ${reads}
        samtools view -@ ${params.threads} -Sb -F 2052 tmp.sam | \
        samtools sort -@ ${params.threads} -o mapped_reads.bam -
        samtools index mapped_reads.bam
      fi

      MEAN_DEPTHS=`samtools depth -aa mapped_reads.bam | awk '{sum+=\$3} END { print int(sum/NR)}'`                                                                                               # We call variants with freyja                                                                                                                                                              if [ \${MEAN_DEPTHS} -gt 20 ]; then
        # We call variants with freyja
 
        freyja variants mapped_reads.bam --minq ${params.freyja_minq} --variants variants_files/test.variants.tsv --depths depth_files/test.depth --ref reference.fasta
        freyja demix variants_files/test.variants.tsv depth_files/test.depth --output demix_files/test.output --confirmedonly --barcodes  /home/external_databases/freyja/RSV_${TYPE}/barcode.csv

        freyja aggregate demix_files/ --output coinfections.tsv
        freyja_status="tak"
        freyja_lineage_1_name=`cat coinfections.tsv  | cut -f3 | tail -1 | cut -d " " -f1 | sed s"|RSVa-||"g | sed s"|RSVb-||"g`
        freyja_lineage_2_name=`cat coinfections.tsv  | cut -f3 | tail -1 | cut -d " " -f2 | sed s"|RSVa-||"g | sed s"|RSVb-||"g` 
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
           freyja_status="blad"
           echo -e "{\\"status\\":\\"\${freyja_status}\\",
                     \\"error_message\\":\\"\${ERR_MSG}\\"}" >> coinfections_freyja.json

         fi
    fi

    """
}
