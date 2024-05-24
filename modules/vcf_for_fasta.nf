process vcf_for_fasta {
    tag "Creating VCF for consensus sequence for sample:\t$sampleId"
    publishDir "${params.results_dir}/${sampleId}/consensus_vcf", mode: 'copy', pattern: "\${vcf_output}*"

    input:
    tuple val(sampleId), path("consensus.fa")
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(sampleId), path("consensus.vcf.gz"), path("consensus.vcf.gz.tbi")

    script:
    """
    N=1 # maksymalna przerwa miedzy mutacjami aby je polaczyc w jeden ciag
    vcf_input="/home/SARS-CoV2/vcf_template/vcf_template.vcf.gz" # sciezka do pliku vcf, paczka vcf z pythona wymaga pliku vcf aby na podstawie tego schematy samemu tworzyc plik output-owy

    # plik ten jest częścią repozytorium w katalogu
    vcf_output="consensus.vcf" # nazwa pliku vcf dla sekwencji konsensusowej
    prep_own_vcf.py ${reference_fasta} consensus.fa \$N \${vcf_input} \${vcf_output}
    bgzip \${vcf_output};
    tabix \${vcf_output}.gz
   """
}
