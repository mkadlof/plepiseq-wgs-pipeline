process vcf_for_fasta {
    tag "vcf_for_fasta:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}/consensus_vcf", mode: 'copy', pattern: "*vcf*"

    input:
    tuple val(sampleId), path("consensus.fa")

    output:
    tuple val(sampleId), path("consensus.vcf.gz"), path("consensus.vcf.gz.tbi")

    script:
    """
    N=1 # maksymalna przerwa miedzy mutacjami aby je polaczyc w jeden ciag
    vcf_input="/home/data/sarscov2/vcf_template/vcf_template.vcf.gz" # sciezka do pliku vcf, paczka vcf z pythona wymaga pliku vcf aby na podstawie tego schematy samemu tworzyc plik output-owy

    # plik ten jest częścią repozytorium w katalogu
    vcf_output="consensus.vcf" # nazwa pliku vcf dla sekwencji konsensusowej
    # usuwamy x ktore teraz w genomie znacza delecje ale ani moj skrypt ani snpeff tego nie zrozumie
    sed -i s'/x//'g consensus.fa
    prep_own_vcf.py \${GENOME_FASTA} consensus.fa \$N \${vcf_input} \${vcf_output}
    bgzip \${vcf_output};
    tabix \${vcf_output}.gz
   """
}
