process modeller {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(target_fasta)

    output:
    tuple val(sampleId), path('alignment.pir'), path('target.B99990001.pdb')

    script:
    """
    if [ -s ${target_fasta} ] ; then
        cp /home/SARS-CoV2/modeller/7dwz.pdb .
        modpy.sh modeller_create_alignment.py ${target_fasta}
        modpy.sh modeller_build_model.py alignment.pir
    else
        echo "Empty fasta file. Skipping modeller."
        touch alignment.pir
        touch target.B99990001.pdb
    fi
    """
}