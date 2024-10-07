process freyja {
    tag "freyja:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "coinfections.tsv"
    containerOptions "--volume ${params.freyja_db_absolute_path_on_host}:/home/external_databases/freyja"

    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), path(ref_genome)

    output:
    tuple val(sampleId), path('coinfections.tsv')

    script:
    """
    count=\$(samtools view "mapped_reads.bam" | wc -l)
    if [ \$count -eq 0 ]; then
        echo "No reads in the bam file"
        touch coinfections.tsv
        echo "No reads mapped to provided reference genome" >> coinfections.tsv
        exit 0
    else
        mkdir variants_files depth_files demix_files
        freyja variants mapped_reads.bam --variants variants_files/test.variants.tsv --depths depth_files/test.depth --ref ${ref_genome}
        freyja demix variants_files/test.variants.tsv depth_files/test.depth --output demix_files/test.output --confirmedonly --barcodes  /home/external_databases/freyja/usher_barcodes.csv
        freyja aggregate demix_files/ --output coinfections.tsv
    fi
    """
}
