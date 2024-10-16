process json_aggregator {
    tag "json_aggregator:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "output.json"

    input:
    val pathogen
    val pipeline_version
    tuple val(sampleId), path(wgsMetrics_json)

    output:
    file('output.json')

    script:
    """
    json_aggregator.py ${pipeline_version} ${pathogen} ${sampleId} \
                       --wgsMetrics ${wgsMetrics_json}
    """
}