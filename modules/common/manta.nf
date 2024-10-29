process introduce_SV_with_manta {
    // This module iterates over all bam_file provided to this module via picard_downsample
    tag "manta:$sampleId"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy'
    container = params.manta_image

    input:
    tuple val(sampleId), path(bam_files), path(bai_files),  path(ref_genome_with_index), path("mediana_per_segment.txt"), val(QC_status_picard), path(consensus_files), val(QC_status_consensus)

    // zarowno zmienna bam_files, bai_files i consnesus_files sa listami
    // moduly wyzej zapewniaja by nazwy plikow byly latwo parsowalne tzn pliki bam maja nazwe downsampling_NAZWA.bam a pliki fasta
    // dla tego segmentu consesnsu_NAZWA.fasta
    // znajac NAZWA mozemy latwo znajdowac odpowiednie pliki

    output:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path(ref_genome_with_index), env(QC_status_exit), emit: fasta_refgenome_and_qc
    tuple val(sampleId), path('dummy_json.json'), emit: json
    // Dal ulatwienia ? na koniec tego segmentu polaczmy wszystkie segmenty w jeden plik, a jelsi dany modul downstream
    // bedzie wymagal sekwencji konkretnego segmentu to tam zrobimy split-a

    script:
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }

    
    """
    ### podiana nazw linkow symbolicznych
    # for link in \$(find . -maxdepth 1 -type l); do
    #    target=\$(readlink "\$link")
    #    base_target=\$(basename "\$target")
    #    mv "\$link" "\$base_target"
    #done

    touch dummy_json.json
    if [[ ${QC_status_picard} == "nie"  && ${QC_status_consensus} == "nie" ]]; then
      
      # both consensus module and picard failed, dummy output
      touch output_consensus_masked_SV.fa
      QC_status_exit="nie"
    elif [[ ${QC_status_picard} == "nie" &&  ${QC_status_consensus} == "tak" ]]; then
      # downsampling failed for ALL segemtns, but consensus produced valid output, ALL input fastas becomes  output files for their respective segemnts
      # We only change the header
      QC_status_exit="tak"
  
      for plik in `ls con*fasta`; do
        plik_new_name=`basename \${plik} ".fasta"`
        HEADER=`head -1 \${plik}`
        cat \${plik} | sed s"|\${HEADER}|\${HEADER}_SV|"g > \${plik_new_name}_SV.fasta

      done
      cat *SV.fasta >> output_consensus_masked_SV.fa
 
    else
      # For, at least, one segment we need to run manta ... however we do not know for which segment
      # for now we analyze each segment separetly (running  manta , introducing mutations and so on)
      # for future reference maybe we can combine these individual steps into one ?

      # We split reference genome into segments, we ensure that after reference_ 
      # segment name is consistence with consensus_ and downsample_ files

      cat ${ref_genome_with_index[final_index]} | awk '{if (substr(\$0, 1, 1)==">") { new_name=\$0; gsub("\\\\.", "_", new_name); gsub("/", "_", new_name);  filename=("reference_"substr(new_name,2) ".fasta"); print \$0 >> filename } else {print toupper(\$0)  >> filename}}'
      
      for bam in `ls *bam`; do
        segment_clean=`basename \${bam}  ".bam" | cut -d "_" -f2- `
        segment_median_coverage=`cat mediana_per_segment.txt | grep -w \${segment_clean} | cut -f2`
        if [ \${segment_median_coverage} -lt 50 ]; then
          # After downsampling segment has poor coverage    
          HEADER=`head -1 consensus_\${segment_clean}.fasta`
          cat consensus_\${segment_clean}.fasta | sed s"|\${HEADER}|\${HEADER}_SV|"g > consensus_\${segment_clean}_SV.fasta
        else
          samtools faidx reference_\${segment_clean}.fasta
          python /opt/docker/manta/bin/configManta.py --bam \${bam} --reference reference_\${segment_clean}.fasta --runDir Manta_results_\${segment_clean}
          python Manta_results_\${segment_clean}/runWorkflow.py -j ${params.threads} --quiet
          if [ -e Manta_results_\${segment_clean}/results/variants/diploidSV.vcf.gz ]; then
            # Manta produced an output for this segment
            # Wywalamy skomplikowane SV jak translokacje itd typun BND
            bcftools view -O z -o manta_results_\${segment_clean}.vcf.gz -i '(FILTER="PASS" | FILTER="MaxDepth" | FILTER="NoPairSupport") && SVTYPE != "BND"' Manta_results_\${segment_clean}/results/variants/diploidSV.vcf.gz
            tabix manta_results_\${segment_clean}.vcf.gz
            ILE_SV=`zcat manta_results_\${segment_clean}.vcf.gz | grep SVTYPE | grep -v INFO | grep -v bcftools | wc -l`

            if [ \${ILE_SV} -gt 0 ]; then
              # MAnta produced at least one  valid SV we must integrate them into input fasta
              cat reference_\${segment_clean}.fasta | bcftools consensus -s - manta_results_\${segment_clean}.vcf.gz  > output_manta_\${segment_clean}.fa
              /home/bin/sarscov2/insert_SV_python2.py consensus_\${segment_clean}.fasta output_manta_\${segment_clean}.fa  consensus_\${segment_clean}_SV.fasta
            else
              HEADER=`head -1 consensus_\${segment_clean}.fasta`
              cat consensus_\${segment_clean}.fasta | sed s"|\${HEADER}|\${HEADER}_SV|"g > consensus_\${segment_clean}_SV.fasta
            fi # koniec if-a na brak SV

          else
            HEADER=`head -1 consensus_\${segment_clean}.fasta`
            cat consensus_\${segment_clean}.fasta | sed s"|\${HEADER}|\${HEADER}_SV|"g > consensus_\${segment_clean}_SV.fasta
          fi # koniec if-a na brak outputu manty
        fi # koniec if-a na zly coverage
      done # koniec petli na teracje po segmentach
      # merge all fasta into a single file
      cat *SV.fasta >> output_consensus_masked_SV.fa  
    fi # koniec if-a na QC
    """
}

