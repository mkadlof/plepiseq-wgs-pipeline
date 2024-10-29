process freyja {
    tag "freyja:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "coinfections.tsv"
    containerOptions "--volume ${params.external_databases_path}:/home/external_databases/"

    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), path(ref_genome_with_index), val(QC_status)

    output:
    tuple val(sampleId), path('coinfections.tsv'), emit: to_pubdir
    tuple val(sampleId), path('coinfections_freyja.json'), emit: json

    script:
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }
    """
    if [ ${QC_status} == "nie" ]; then
      # echo "No reads in the bam file"
      touch coinfections.tsv
      touch coinfections_freyja.json
    else
      mkdir variants_files depth_files demix_files
      freyja variants mapped_reads.bam --variants variants_files/test.variants.tsv --depths depth_files/test.depth --ref ${ref_genome_with_index[final_index]}
      freyja demix variants_files/test.variants.tsv depth_files/test.depth --output demix_files/test.output --confirmedonly --barcodes  /home/external_databases/freyja/usher_barcodes.csv
      freyja aggregate demix_files/ --output coinfections.tsv
      touch coinfections_freyja.json
    fi
    """
}
