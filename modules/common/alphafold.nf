process alphafold {
    tag "alphafold:${sampleId}"
    maxRetries 3
    cpus 15 
    // Set to 15. We have 96 COU on compute and  8 GPUs,  so 15 CPUs per task  * 7 tasks will require 95 CPus total 
    errorStrategy 'retry' 
    maxForks 8 
    // Let us leave one GPU free just in case a "failed" process rebounce and quickly ask for an empy GPU
    // But for some reason nextflow executes not "maxforks" but "maxforks -1" processes
    container  = params.alphafold_image
    containerOptions "--volume ${params.alphafold_databases_path}:/db --gpus all"
    // publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}_spike.pdb"

    input:
    tuple val(sampleId), path("*fasta"), val(QC_status)

    output:
    tuple val(sampleId), path("${sampleId}_spike.pdb"), emit: to_pubdir
    // tuple val(sampleId), path('modeller.json'), emit: json

    script:
    """
    for link in \$(find . -maxdepth 1 -type l); do
        target=\$(readlink "\$link")
        base_target=\$(basename "\$target")
        mv "\$link" "\$base_target"
    done

    if [[ ${QC_status} == "nie" || ${params.species} != "SARS-CoV-2" ]]; then
      touch ${sampleId}_spike.pdb
      echo -e "{\\"protein_structure_status\\":\\"nie\\"}" >> modeller.json
    else
      if [ -e "nextalign_gene_S.translation.fasta" ]; then
        # There should be a limit on number of "-" and "X" in a sequence 
       
        cat nextalign_gene_S.translation.fasta | tr -d "-" | tr -d "X" | tr -d "*" >> tmp
        mv tmp nextalign_gene_S.translation.fasta
        target_fasta="nextalign_gene_S.translation.fasta"

        sed -i s"|n_cpu: int = 8|n_cpu: int = ${task.cpus}|"g /app/alphafold/alphafold/data/tools/jackhmmer.py
        sed -i s"|n_cpu: int = 4|n_cpu: int = ${task.cpus}|"g /app/alphafold/alphafold/data/tools/hhblits.py

        # update config to run only model_1 from alphafold config
        sed -i -E "s|'model_1',|['model_1']|g" /app/alphafold/alphafold/model/config.py
        sed -i -zE "s|'model_2',\\s*'model_3',\\s*'model_4',\\s*'model_5',\\s*||g" /app/alphafold/alphafold/model/config.py
        # remove models 2 .. .5 from config

        sleep `python -c 'import random; print(random.randint(1, 15))'`
        nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits > tmp_smi.log

        # Determine which GPU is free 
        ID=0
        while read L; do
          SMI_1=`echo \${L} | cut -d " " -f1`

          if [ \${SMI_1} -lt 50 ]; then
            break
          else
            ID=`echo "\${ID} + 1" | bc -l`
          fi
       done < tmp_smi.log

       if [ \${ID} -gt 7 ]; then
         exit 1
       fi
       # Steer alphafold to use free gpu
       export CUDA_VISIBLE_DEVICES=\${ID}
 
       mkdir wyniki
       python /app/alphafold/run_alphafold.py  --fasta_paths="\${target_fasta}" \
                                               --data_dir="/db/" \
                                               --db_preset="reduced_dbs" \
                                               --output_dir="wyniki" \
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
        # --uniref30_database_path="/db/uniref30/UniRef30_2021_03" \

        sleep `python -c 'import random; print(random.randint(5, 8))'`
        cp wyniki/`basename \${target_fasta} ".fasta"`/ranked_0.pdb ${sampleId}_spike.pdb
        
        # align all proteins to a common reference
        # TO DO
      else
        touch ${sampleId}_spike.pdb
      fi
    fi
   """
}
