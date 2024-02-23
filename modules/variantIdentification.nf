process variantIdentification {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(freebayes_masked_fa), path(lofreq_masked_fa), path(varscan_masked_fa), path(consensus_masked_fa)

    output:
    tuple val(sampleId), path("nextclade_linages"), path('nextstrain_lineage.csv'), path('pangolin_lineage.csv')

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
                  --output-all nextclade_linages \
                  all_sequences.fa
    """
}