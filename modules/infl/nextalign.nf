process nextalign {
    tag "nextalign:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "consensus_M2.fasta"

    input:
    tuple val(sampleId), env("REF_GENOME_ID_MINI")
    tuple val(sampleId), path("consensus_*.fasta")
    val("NEXTALIGN_DB")

    output:
    tuple val(sampleId), path("nextalign_gene_*.translation.fasta")
    tuple val(sampleId), path("consensus_M2.fasta")

    script:
    """
    # === This workaround because nextflow replace glob patterns with numbers ===
    # It goal is to rename symlinks to their targets names.
    for link in \$(find . -maxdepth 1 -type l); do
        target=\$(readlink "\$link")
        base_target=\$(basename "\$target")
        mv "\$link" "\$base_target"
    done
    # === End of workaround ===

    echo "REF_GENOME_ID_MINI: \${REF_GENOME_ID_MINI}"

    TO_NEXTALIGN='H3N2 H5N1 H5N5 H5N6 H5N8 H7N9'
    if [[ \${TO_NEXTALIGN[@]} =~ \${REF_GENOME_ID_MINI} ]]; then
        for FILE in `ls *.fasta`
        do
            SEGMENT=`echo \${FILE} | awk -F "[._]" '{print \$2}'`
            nextalign run -r ${NEXTALIGN_DB}/\${REF_GENOME_ID_MINI}/\${SEGMENT}.fasta \
                          -m ${NEXTALIGN_DB}/\${REF_GENOME_ID_MINI}/\${REF_GENOME_ID_MINI}.gff \
                          -g \${SEGMENT} \
                          -O . \${FILE}
            if [[ \${SEGMENT} == 'NA' || \${SEGMENT} == 'PA' ]]; then
                cat nextalign_gene_\${SEGMENT}.translation.fasta >> to_resistance.fasta
            fi

            if [ \${SEGMENT} == 'MP' ]; then
                prep_M2.py \${REF_GENOME_ID_MINI} MP.fasta ${NEXTALIGN_DB} . .
                nextalign run -r M2.fasta \
                              -m ${NEXTALIGN_DB}/\${REF_GENOME_ID_MINI}/\${REF_GENOME_ID_MINI}.gff \
                              -g M2 \
                              -O . consensus_M2.fasta
            fi
        done
    fi
    """
}
