process json_aggregator_sars_illumina {
    tag "json_aggregator:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}.json"
    container  = params.main_image
    containerOptions "--volume ${params.projectDir}:/home/projectDir:ro"

    input:
    tuple val(sampleId), path(fastqc_pre_json_forward), path(fastqc_pre_json_reverse),
          path(kraken_contamination),
          path(fastqc_post_json_forward), path(fastqc_post_json_reverse),
          path(freyja),
          path(coinfection),
          path(dehumanized),
          path(wgsMetrics),
          path(consensus_json),
          path(pangolin_json), 
          path(nextclade_json),
          path(snpeff),
          path(modeller)


    output:
    path("${sampleId}.json")

    script:
    """
    branch=master
    version=\$(cat /home/projectDir/.git/refs/heads/\${branch})
    version=\${version:0:7}
 
    json_aggregator.py  --version \${version} \
                        --pathogen "${params.species}" \
                        --sampleId "${sampleId}" \
                        --fastqc_pre "${fastqc_pre_json_forward}" "${fastqc_pre_json_reverse}" \
                        --fastqc_post "${fastqc_post_json_forward}" "${fastqc_post_json_reverse}" \
                        --contamination "${kraken_contamination}" \
                        --freyja "${freyja}" \
                        --coinfection "${coinfection}" \
                        --dehumanized "${dehumanized}" \
                        --wgsMetrics "${wgsMetrics}"  \
                        --consensus "${consensus_json}" \
                        --pangolin "${pangolin_json}" \
                        --nextclade "${nextclade_json}" \
                        --modeller "${modeller}" \
                        --snpeff ${snpeff}

    cp output.json ${sampleId}.json
    """
}

process json_aggregator_rsv_illumina {
    tag "json_aggregator:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}.json"
    container  = params.main_image
    containerOptions "--volume ${params.projectDir}:/home/projectDir:ro"

    input:
    tuple val(sampleId), path(fastqc_pre_json_forward), path(fastqc_pre_json_reverse),
          path(kraken_contamination),
          path(fastqc_post_json_forward), path(fastqc_post_json_reverse),
          path(dehumanized),
          path(wgsMetrics),
          path(consensus_json),
          path(pangolin_json),
          path(nextclade_json),
          path(snpeff),
          path(modeller)


    output:
    path("${sampleId}.json")

    script:
    """
    branch=master
    version=\$(cat /home/projectDir/.git/refs/heads/\${branch})
    version=\${version:0:7}
    json_aggregator.py  --version \${version} \
                        --pathogen "${params.species}" \
                        --sampleId "${sampleId}" \
                        --fastqc_pre "${fastqc_pre_json_forward}" "${fastqc_pre_json_reverse}" \
                        --fastqc_post "${fastqc_post_json_forward}" "${fastqc_post_json_reverse}" \
                        --contamination "${kraken_contamination}" \
                        --dehumanized "${dehumanized}" \
                        --wgsMetrics "${wgsMetrics}"  \
                        --consensus "${consensus_json}" \
                        --pangolin "${pangolin_json}" \
                        --nextclade "${nextclade_json}" \
                        --snpeff ${snpeff}

    cp output.json ${sampleId}.json
    """
}


process json_aggregator_influenza_illumina {
    tag "json_aggregator:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}.json"
    container  = params.main_image
    containerOptions "--volume ${params.projectDir}:/home/projectDir:ro"

    input:
    tuple val(sampleId), path(fastqc_pre_json_forward), path(fastqc_pre_json_reverse),
          path(kraken_contamination),
          path(fastqc_post_json_forward), path(fastqc_post_json_reverse),
          path(reassortment_json),
          path(dehumanized),
          path(wgsMetrics),
          path(consensus_json),
          path(pangolin_json),
          path(nextclade_json),
          path(snpeff),
          path(modeller),
          path(resistance_json)


    output:
    path("${sampleId}.json")

    script:
    """
    branch=master
    version=\$(cat /home/projectDir/.git/refs/heads/\${branch})
    version=\${version:0:7}

    touch ${sampleId}.json
    """
}


process json_aggregator_sars_nanopore {
    tag "json_aggregator:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}.json"
    container  = params.main_image
    containerOptions "--volume ${params.projectDir}:/home/projectDir:ro"

    input:
    tuple val(sampleId), path(fastqc_pre_json_forward),
          path(kraken_contamination),
          path(freyja),
          path(coinfection),
          path(dehumanized),
          path(wgsMetrics),
          path(consensus_json),
          path(pangolin_json),
          path(nextclade_json),
          path(snpeff),
          path(modeller)


    output:
    path("${sampleId}.json")

    script:
    """
    branch=master
    version=\$(cat /home/projectDir/.git/refs/heads/\${branch})
    version=\${version:0:7}
    
    json_aggregator.py  --version \${version} \
                        --pathogen "${params.species}" \
                        --sampleId "${sampleId}" \
                        --fastqc_pre "${fastqc_pre_json_forward}" \
                        --contamination "${kraken_contamination}" \
                        --freyja "${freyja}" \
                        --coinfection "${coinfection}" \
                        --dehumanized "${dehumanized}" \
                        --wgsMetrics "${wgsMetrics}"  \
                        --consensus "${consensus_json}" \
                        --pangolin "${pangolin_json}" \
                        --nextclade "${nextclade_json}" \
                        --modeller "${modeller}" \
                        --snpeff ${snpeff}

    cp output.json ${sampleId}.json
    """
}

process json_aggregator_rsv_nanopore {
    tag "json_aggregator:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}.json"
    container  = params.main_image
    containerOptions "--volume ${params.projectDir}:/home/projectDir:ro"

    input:
    tuple val(sampleId), path(fastqc_pre_json_forward),
          path(kraken_contamination),
          path(dehumanized),
          path(wgsMetrics),
          path(consensus_json),
          path(pangolin_json),
          path(nextclade_json),
          path(snpeff),
          path(modeller)

    output:
    path("${sampleId}.json")

    script:
    """
    branch=master
    version=\$(cat /home/projectDir/.git/refs/heads/\${branch})
    version=\${version:0:7}
    json_aggregator.py  --version \${version} \
                        --pathogen "${params.species}" \
                        --sampleId "${sampleId}" \
                        --fastqc_pre "${fastqc_pre_json_forward}" \
                        --contamination "${kraken_contamination}" \
                        --dehumanized "${dehumanized}" \
                        --wgsMetrics "${wgsMetrics}"  \
                        --consensus "${consensus_json}" \
                        --pangolin "${pangolin_json}" \
                        --nextclade "${nextclade_json}" \
                        --snpeff ${snpeff}

    cp output.json ${sampleId}.json
    
    """
}

process json_aggregator_influenza_nanopore {
    tag "json_aggregator:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}.json"
    container  = params.main_image
    containerOptions "--volume ${params.projectDir}:/home/projectDir:ro"

    input:
    tuple val(sampleId), path(fastqc_pre_json_forward),
          path(kraken_contamination),
          path(reassortment_json),
          path(dehumanized),
          path(wgsMetrics),
          path(consensus_json),
          path(pangolin_json),
          path(nextclade_json),
          path(snpeff),
          path(modeller),
          path(resistance_json)

    output:
    path("${sampleId}.json")

    script:
    """
    branch=master
    version=\$(cat /home/projectDir/.git/refs/heads/\${branch})
    version=\${version:0:7}

    touch ${sampleId}.json
    """
}

