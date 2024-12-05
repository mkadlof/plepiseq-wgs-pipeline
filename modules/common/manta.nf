process introduce_SV_with_manta {
    // This module iterates over all bam_file provided to this module via picard_downsample
    tag "manta:$sampleId"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "consensus*fasta"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "output_consensus_masked_SV.fa"
    container = params.manta_image

    input:
    tuple val(sampleId), path(bam_files), path(bai_files),  path(ref_genome_with_index), path("mediana_per_segment.txt"), val(QC_status_picard), path(consensus_files), val(QC_status_consensus)

    // zarowno zmienna bam_files, bai_files i consnesus_files sa listami
    // moduly wyzej zapewniaja by nazwy plikow byly latwo parsowalne tzn pliki bam maja nazwe downsampling_NAZWA.bam a pliki fasta
    // dla tego segmentu consesnsu_NAZWA.fasta
    // znajac NAZWA mozemy latwo znajdowac odpowiednie pliki

    output:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path(ref_genome_with_index), env(QC_status_exit), emit: fasta_refgenome_and_qc
    tuple val(sampleId), path('consensus.json'), emit: json
    tuple val(sampleId), path('consensus_*.fasta'), emit: to_pubdir
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
    export LC_ALL=en_US.utf-8
    if [[ ${QC_status_picard} == "nie"  && ${QC_status_consensus} == "nie" ]]; then
      
      # both consensus module and picard failed, dummy output
      touch output_consensus_masked_SV.fa
      touch consensus_dummy.fasta
      QC_status_exit="nie"
      ERR_MSG="Failed QC"
      python3 /home/parse_make_consensus.py --status "nie" --error "\${ERR_MSG}" -o consensus.json

    elif [[ ${QC_status_picard} == "nie" &&  ${QC_status_consensus} == "tak" ]]; then
      # downsampling failed for ALL segemtns, but consensus produced valid output, ALL input fastas becomes  output files for their respective segemnts
      # We only change the header
      QC_status_exit="tak"
  
      for plik in `ls con*fasta`; do
        plik_new_name=`basename \${plik} ".fasta"`
        HEADER=`head -1 \${plik}`
        cat \${plik} | sed s"|\${HEADER}|\${HEADER}_SV|"g > \${plik_new_name}_SV.fasta

      done
      
      # make json
      ls *_SV.fasta | tr " " "\\n" >> list_of_fasta.txt
      python3 /home/parse_make_consensus.py --status "tak" -o consensus.json --input_fastas list_of_fasta.txt --output_path "${params.results_dir}/${sampleId}"
      
      # merge fastas of individual segments to a single file
      cat *SV.fasta >> output_consensus_masked_SV.fa

    else
      QC_status_exit="tak"
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
            # Wwywalamy rowniez takie SV ktore nie sa homozygotyczne, niepwene mutacje z GT 0/1 beda usuwane 
            bcftools view -O z -o manta_results_\${segment_clean}.vcf.gz -i '(FILTER="PASS" | FILTER="MaxDepth" | FILTER="NoPairSupport") && SVTYPE != "BND" && GT!="het"' Manta_results_\${segment_clean}/results/variants/diploidSV.vcf.gz
            tabix manta_results_\${segment_clean}.vcf.gz
            ILE_SV=`zcat manta_results_\${segment_clean}.vcf.gz | grep SVTYPE | grep -v INFO | grep -v bcftools | wc -l`

            if [ \${ILE_SV} -gt 0 ]; then
              # MAnta produced at least one  valid SV we must integrate them into input fasta
              cat reference_\${segment_clean}.fasta | bcftools consensus -s - manta_results_\${segment_clean}.vcf.gz  > output_manta_\${segment_clean}.fa
              # add _manta do fasta header for manta
              HEADER=`head -1 output_manta_\${segment_clean}.fa`
              sed -i s"|\${HEADER}|\${HEADER}_manta|"g output_manta_\${segment_clean}.fa
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

      # create json
      ls *_SV.fasta | tr " " "\\n" >> list_of_fasta.txt
      python3 /home/parse_make_consensus.py --status "tak" -o consensus.json --input_fastas list_of_fasta.txt --output_path "${params.results_dir}/${sampleId}"
      
      # merge all fasta into a single file
      cat *SV.fasta >> output_consensus_masked_SV.fa 
     fi # koniec if-a na QC

    # in consensus json remove _SV from segment name
    sed -i s'|_SV"|"|'g consensus.json 
    # Fasta headears should follow specific naming scheme "{segment_name}|{sample_name}" but this BREAKS downstream modules
    # Thus to PUBDIR we will push different files than the one pushed to downstream modules 
    # TO DO
    """
}

