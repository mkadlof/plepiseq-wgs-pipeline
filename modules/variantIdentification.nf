process variantIdentification {
    publishDir "results/${sampleId}", mode: 'symlink'
    containerOptions "--volume ${params.pangolin_db_absolute_path_on_host}:/home/SARS-CoV2/pangolin"

    input:
    tuple val(sampleId), path(freebayes_masked_fa), path(lofreq_masked_fa), path(varscan_masked_fa), path(consensus_masked_fa)

    output:
    tuple val(sampleId), path("nextclade_lineages"), path('nextstrain_lineage.csv'), path('pangolin_lineage.csv')

    script:
    """
    cat ${consensus_masked_fa} \
        ${freebayes_masked_fa} \
        ${lofreq_masked_fa} \
        ${varscan_masked_fa} > all_sequences.fa
    pangolin --outfile pangolin_lineage.csv \
             --threads ${params.threads} \
             all_sequences.fa
    nextclade run --input-dataset ${params.nextclade_db_dir} \
                  --output-csv nextstrain_lineage.csv \
                  --output-all nextclade_lineages \
                  all_sequences.fa
    """
}