process modeller {
    tag "modeller:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}_spike.pdb"

    input:
    tuple val(sampleId), path(target_fasta)

    output:
    tuple val(sampleId), path('alignment.pir'), path("${sampleId}_spike.pdb")

    script:
    """
    if [ -s ${target_fasta} ] ; then
        cp /home/SARS-CoV2/modeller/7dwz.pdb .
        modpy.sh modeller_create_alignment.py ${target_fasta}
        modpy.sh modeller_build_model.py alignment.pir
        cp target.B99990001.pdb ${sampleId}_spike.pdb 
    else
        echo "Empty fasta file. Skipping modeller."
        touch alignment.pir
        touch ${sampleId}_spike.pdb
    fi
    """
}
