process make_genome_from_vcf {
    tag "varScan:${sampleId}"
    container  = params.main_image

    input:
    tuple val(sampleId),  path('genome.fasta'),  val(QC_status), path('input.vcf.gz'), path('input.vcf.gz.tbi')

    output:
    tuple val(sampleId), path('sample_genome.fa'), val(QC_status), emit: fasta

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch varscan.fa
      touch detected_variants_varscan.txt
    else
  
      cat genome.fasta | bcftools consensus input.vcf.gz > sample_genome.fa
    fi
    """
}
