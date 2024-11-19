process json_aggregator {
    tag "json_aggregator:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "output.json"
    container  = params.main_image
    containerOptions "--volume ${params.projectDir}:/home/projectDir:ro"

    input:
    tuple val(sampleId), path(wgsMetrics_json),
        path(segment_bedgraphs_files_txt),
        path(consensus_json), path(list_of_fasta_files),
        path(kraken_contamination),
        path(fastqc_pre_json_forward), path(fastqc_pre_json_reverse),
        path(fastqc_post_json_forward), path(fastqc_post_json_reverse),
        path(pangolin_json), path(nextclade_json)

    output:
    file('output.json')

    script:
    """
    branch=master
    version=\$(cat /home/projectDir/.git/refs/heads/\${branch})
    version=\${version:0:7}
    json_aggregator.py \${version} ${params.species} ${sampleId} ${params.results_dir}/${sampleId} \
                       --wgsMetrics ${wgsMetrics_json} \
                       --segment_bedgraphs_files ${segment_bedgraphs_files_txt} \
                       --consensus ${consensus_json} \
                       --list_of_fasta_files ${list_of_fasta_files} \
                       --contamination ${kraken_contamination} \
                       --fastqc_pre ${fastqc_pre_json_forward} ${fastqc_pre_json_reverse} \
                       --fastqc_post ${fastqc_post_json_forward} ${fastqc_post_json_reverse} \
                       ${!pangolin_json.contains('/non-existent') ? "--pangolin ${pangolin_json}" : ""} \
                       --nextclade ${nextclade_json}
    """
}

process json_aggregator_nanopore {
    tag "json_aggregator:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "output.json"
    container  = params.main_image
    containerOptions "--volume ${params.projectDir}:/home/projectDir:ro"

    input:
    tuple val(sampleId), path(wgsMetrics_json),
        path(segment_bedgraphs_files_txt),
        path(consensus_json), path(list_of_fasta_files),
        path(kraken_contamination),
        path(fastqc_pre_json_forward),
        path(pangolin_json), path(nextclade_json)

    output:
    file('output.json')

    script:
    """
    branch=master
    version=\$(cat /home/projectDir/.git/refs/heads/\${branch})
    version=\${version:0:7}
    json_aggregator.py \${version} ${params.species} ${sampleId} ${params.results_dir}/${sampleId} \
                       --wgsMetrics ${wgsMetrics_json} \
                       --segment_bedgraphs_files ${segment_bedgraphs_files_txt} \
                       --consensus ${consensus_json} \
                       --list_of_fasta_files ${list_of_fasta_files} \
                       --contamination ${kraken_contamination} \
                       --fastqc_pre ${fastqc_pre_json_forward} \
                       --pangolin ${pangolin_json} \
                       --nextclade ${nextclade_json}
    """
}

