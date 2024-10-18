process json_aggregator {
    tag "json_aggregator:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "output.json"

    input:
    val pathogen
    val pipeline_version
    tuple val(sampleId), path(wgsMetrics_json), path(segment_bedgraphs_files_txt)

    output:
    file('output.json')

    script:
    """
    json_aggregator.py ${pipeline_version} ${pathogen} ${sampleId} ${params.results_dir}/${sampleId}\
                       --wgsMetrics ${wgsMetrics_json} \
                       --segment_bedgraphs_files ${segment_bedgraphs_files_txt}
    """
}