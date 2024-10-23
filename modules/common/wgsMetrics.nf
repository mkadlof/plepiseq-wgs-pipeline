process picard_wgsMetrics {
    tag "wgsMetrics:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "picard_statistics.txt"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "segment_*.bedgraph"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "coverage_histogram_data.tsv"

    input:
    tuple val(sampleId), path(bam), path(bai), val(QC_status), path(ref_genome)

    output:
    tuple val(sampleId), path('picard_statistics.txt'), emit: statistics
    tuple val(sampleId), path("wgsMetrics.json"), path("segment_bedgraphs_files.txt"), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch picard_statistics.txt
      touch wgsMetrics.json
      touch segment_bedgraphs_files.txt
    else
      java -jar /opt/picard/picard.jar CollectWgsMetrics --REFERENCE_SEQUENCE ${ref_genome} \
                                                         --MINIMUM_BASE_QUALITY ${params.quality_initial} \
                                                         --MINIMUM_MAPPING_QUALITY ${params.min_mapq} \
                                                         --INPUT ${bam} \
                                                         --OUTPUT picard_statistics.txt


      bedtools genomecov -d -ibam ${bam} > genomecov.bedgraph

      cat  genomecov.bedgraph | awk '{
      gsub("/", "_"); print
      }' >> genomecov_fix.bedgraph

      # Split the genomecov.bedgraph file into segments using awk
      input_file="genomecov_fix.bedgraph"
      output_prefix="segment_"
      awk '{
          if (\$1 != prev) {
              if (prev != "") {
                  close(out)
              }
              prev = \$1
              out = "'"\$output_prefix"'" prev ".bedgraph"
          }
          print >> out
      }' "\$input_file"

      rm genomecov.bedgraph genomecov_fix.bedgraph
      ls segment_*.bedgraph > segment_bedgraphs_files.txt

      parse_wgsMetrics.py picard_statistics.txt wgsMetrics.json
    fi
    """
}
