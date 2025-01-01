process alphafold {
    tag "alphafold:${sampleId}"
    maxRetries 3
    cpus 15 
    // Set to 15. We have 96 COU on compute and  8 GPUs,  so 15 CPUs per task  * 7 tasks will require 95 CPus total 
    errorStrategy 'retry' 
    maxForks 8 
    // Let us leave one GPU free just in case a "failed" process rebounce and quickly ask for an empy GPU
    container  = params.alphafold_image
    containerOptions "--volume ${params.alphafold_databases_path}:/db --gpus all"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "*.pdb"

    input:
    tuple val(sampleId), path("*fasta"), val(QC_status)

    output:
    tuple val(sampleId), path("*.pdb"), emit: to_pubdir // return any number of pdbs produce by this module
    tuple val(sampleId), path('alphafold.json'), emit: json

    script:
    """

    run_alpfafold() {
      # Zmienna $1 to plik z fasta
      # Zmienna $2 to sciezka do katalogu z wynikami
      python /app/alphafold/run_alphafold.py  --fasta_paths="\$1" \
                                               --data_dir="/db/" \
                                               --db_preset="reduced_dbs" \
                                               --output_dir="\$2" \
                                               --uniref90_database_path="/db/uniref50/uniref50.fasta" \
                                               --mgnify_database_path="/db/mgnify/mgy_clusters_2022_05.fa" \
                                               --small_bfd_database_path="/db/small_bfd/bfd-first_non_consensus_sequences.fasta" \
                                               --template_mmcif_dir="/db/pdb_mmcif/mmcif_files/" \
                                               --max_template_date="2024-05-14" \
                                               --obsolete_pdbs_path="/db/pdb_mmcif/obsolete.dat" \
                                               --use_gpu_relax=true \
                                               --pdb70_database_path="/db/pdb70/pdb70" \
                                               --models_to_relax=best \
                                               --model_preset=monomer

     # Options that must be set when  db_preset is set to "full_dbs" , which is a default setting

        # --bfd_database_path="/db/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
        # --uniref30_database_path="/db/uniref30/UniRef30_2021_03" 
    }
   
    # Restore names of fasta files as produced by nextalign module
    for link in \$(find . -maxdepth 1 -type l); do
        target=\$(readlink "\$link")
        base_target=\$(basename "\$target")
        mv "\$link" "\$base_target"
    done

    # in case multiple processes are spawned by nextflow at the same time add random sleep...
    sleep `python -c 'import random; print(random.randint(10, 60))'`

    # increase number of CPUs for jackhammer and hhblits for alignment
    sed -i s"|n_cpu: int = 8|n_cpu: int = ${task.cpus}|"g /app/alphafold/alphafold/data/tools/jackhmmer.py
    sed -i s"|n_cpu: int = 4|n_cpu: int = ${task.cpus}|"g /app/alphafold/alphafold/data/tools/hhblits.py

    # update alphafold config to run only "model_1"
    sed -i -E "s|'model_1',|['model_1']|g" /app/alphafold/alphafold/model/config.py
    sed -i -zE "s|'model_2',\\s*'model_3',\\s*'model_4',\\s*'model_5',\\s*||g" /app/alphafold/alphafold/model/config.py

    # Determine which GPU is free
    nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits > tmp_smi.log
    ID=0
    while read L; do
      SMI_1=`echo \${L} | cut -d " " -f1`
      if [ \${SMI_1} -lt 50 ]; then
        break
      else
        ID=`echo "\${ID} + 1" | bc -l`
      fi
    done < tmp_smi.log

    # If no device is free than send exit code 1 and retry
    if [ \${ID} -gt 7 ]; then
       exit 1
    fi
    
    ## Steer alphafold to use available gpu
    export CUDA_VISIBLE_DEVICES=\${ID}
  
    # directory for alphafold output
    mkdir wynik

    if [ ${QC_status} == "nie" ]; then
      # failed QC 
      touch ${sampleId}.pdb
      ERR_MSG="This modue was entered with bad QC"
      echo -e "{\\"status\\":\\"nie\\", 
                \\"error_message\\": \\"\${ERR_MSG}\\"}" >> alphafold.json
    elif [ ${params.species} == "RSV" ]; then
      # For RSV we predict structures of G and F proteins

      cat nextalign_gene_F.translation.fasta | tr -d "-" | tr -d "X" | tr -d "*" >> tmp
      mv tmp nextalign_gene_F.translation.fasta
      target_fasta_F="nextalign_gene_F.translation.fasta"
      run_alpfafold "\${target_fasta_F}" wynik
      
      cp wynik/`basename \${target_fasta_F} ".fasta"`/ranked_0.pdb ${sampleId}_F.pdb
      # align all proteins to a common reference
      # TO DO

      # json
      pdb_path_1="${params.results_dir}/${sampleId}/${sampleId}_F.pdb"


      cat nextalign_gene_G.translation.fasta | tr -d "-" | tr -d "X" | tr -d "*" >> tmp
      mv tmp nextalign_gene_G.translation.fasta
      target_fasta_G="nextalign_gene_G.translation.fasta"
      run_alpfafold "\${target_fasta_G}" wynik
      cp wynik/`basename \${target_fasta_G} ".fasta"`/ranked_0.pdb ${sampleId}_G.pdb
      # againg align to reference
      
       
      pdb_path_2="${params.results_dir}/${sampleId}/${sampleId}_G.pdb"
      echo -e "{\\"protein_structure_status\\":\\"tak\\",
                \\"protein_structure_data\\":[{\\"protein_name\\":\\"F\\",
                                               \\"pdb_file\\":\\"\${pdb_path_1}\\"
                                              },
                                              {
                                              \\"protein_name\\":\\"G\\",
                                               \\"pdb_file\\":\\"\${pdb_path_2}\\"
                                              }
                                              ]}" >> alphafold.json

    elif [ ${params.species} == "Influenza" ]; then
      cat nextalign_gene_HA.translation.fasta | tr -d "-" | tr -d "X" | tr -d "*" >> tmp
      mv tmp nextalign_gene_HA.translation.fasta
      target_fasta_HA="nextalign_gene_HA.translation.fasta"
      run_alpfafold "\${target_fasta_HA}" wynik

      cp wynik/`basename \${target_fasta_HA} ".fasta"`/ranked_0.pdb ${sampleId}_HA.pdb

      # align all proteins to a common reference
      # TO DO
      # json
      pdb_path_1="${params.results_dir}/${sampleId}/${sampleId}_HA.pdb"
     
      cat nextalign_gene_NA.translation.fasta | tr -d "-" | tr -d "X" | tr -d "*" >> tmp
      mv tmp nextalign_gene_NA.translation.fasta
      target_fasta_NA="nextalign_gene_NA.translation.fasta"
      run_alpfafold "\${target_fasta_NA}" wynik

      cp wynik/`basename \${target_fasta_NA} ".fasta"`/ranked_0.pdb ${sampleId}_NA.pdb

      # align all proteins to a common reference
      # TO DO
      # json
      pdb_path_2="${params.results_dir}/${sampleId}/${sampleId}_HA.pdb"
      echo -e "{\\"protein_structure_status\\":\\"tak\\",
                \\"protein_structure_data\\":[{\\"protein_name\\":\\"HA\\",
                                               \\"pdb_file\\":\\"\${pdb_path_1}\\"
                                              },
                                              {
                                              \\"protein_name\\":\\"NA\\",
                                               \\"pdb_file\\":\\"\${pdb_path_2}\\"
                                              }
                                              ]}" >> alphafold.json

 
    elif [ ${params.species} == "SARS-CoV-2" ]; then
      if [ -e "nextalign_gene_S.translation.fasta" ]; then
       
        cat nextalign_gene_S.translation.fasta | tr -d "-" | tr -d "X" | tr -d "*" >> tmp
        mv tmp nextalign_gene_S.translation.fasta
        target_fasta="nextalign_gene_S.translation.fasta" 
        run_alpfafold "\${target_fasta}" wynik

        # Give some time to clear up memory from a device ...
        sleep `python -c 'import random; print(random.randint(4, 35))'`
        cp wynik/`basename \${target_fasta} ".fasta"`/ranked_0.pdb ${sampleId}_spike.pdb
        
        # align all proteins to a common reference
        # TO DO

        # json
        pdb_path="${params.results_dir}/${sampleId}/${sampleId}_spike.pdb"
        echo -e "{\\"protein_structure_status\\":\\"tak\\",
                \\"protein_structure_data\\":[{\\"protein_name\\":\\"Spike\\",
                                               \\"pdb_file\\":\\"\${pdb_path}\\"
                                               }]}" >> alphafold.json
      else
        touch ${sampleId}.pdb
        ERR_MSG="No fasta for Spike protein"
        echo -e "{\\"status\\":\\"nie\\",
                  \\"error_message\\": \\"\${ERR_MSG}\\"}" >> alphafold.json

      fi
    fi
   """
}
