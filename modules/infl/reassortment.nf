process reassortment {
    tag "reassortment:${sampleId}"
    maxForks 5

    input:
    tuple val(sampleId), path('subtype_scores_each_segment.txt'), path('subtype_counts_each_segment.txt'), env("REF_GENOME_ID")

    output:
    tuple val(sampleId), path('hybrid_genome.fasta')
    tuple val(sampleId), path('hybrid_primers.bed')

    script:
    """
    find_segement_position() {
        local VALUES=("\${!1}")
        local hit=\$2
        local i=1
        for ela in \${VALUES[@]}
        do
            if [ "\$ela" == "\$hit" ]; then
                echo \${i}
                return 0
            else
                ((i++))
            fi
        done
    }

    ALL_GENOMES=(`ls /home/data/infl/genomes`)
    ALL_SEGMENTS=(PB2 PB1 PA HA NP NA MP NS)

    cat /home/data/infl/genomes/\${REF_GENOME_ID}/\${REF_GENOME_ID}.fasta | \
        awk -v ID=\${REF_GENOME_ID} '{if (substr(\$0, 1, 1)==">") {filename=(ID"_"substr(\$0,2) ".fasta"); print \$0"_"ID >> filename } else {print toupper(\$0)  >> filename}}'

    # For each segment, we analyze the file subtype_scores_each_segment.txt with a matrix where the
    # rows are the names of subtypes, the columns are the segments.
    # These 5 lists should hold, in order: the segment identifier, the subtype with the highest
    # alignment score, the score of the expected subtype / score of the best subtype, the
    # similarity between the expected subtype and the found subtype, and the average coverage on
    # the best segment.

    SEGMENTS_REASSORTMENT=()
    FOUND_SUBTYPES_REASSORTMENT=()
    FOUND_SUBTYPES_SCORE_RATIO=()
    FOUND_SUBTYPES_SEQ_SIMILARITY=()
    FOUND_SUBTYPES_COUNTS=()

    for segment in \${ALL_SEGMENTS[@]}; do
        SEGMENTS_REASSORTMENT+=(\${segment})

        # We need to extract data from a file so first we determine in which column the data
        # for a given segment is located.
        LIST_OF_IDS=`cat subtype_scores_each_segment.txt | head -1`
        SEGMENT_POS=`find_segement_position LIST_OF_IDS[@] \${segment}`

        LIST_OF_IDS_COUNTS=`cat subtype_counts_each_segment.txt | head -1`
        SEGMENT_POS_COUNTS=`find_segement_position LIST_OF_IDS_COUNTS[@] \${segment}`

        # name of the subtype with the best result for this segment
        SEGMENT_best=`cat subtype_scores_each_segment.txt | sort -rnk\${SEGMENT_POS} | head -1 | cut -f1`

        # what is the alignment score for this segment from this subtype
        SEGMENT_best_score=`cat subtype_scores_each_segment.txt | sort -rnk\${SEGMENT_POS} | head -1 | cut -f\${SEGMENT_POS}`

        # what is the score for this segment of the sequence found based on the mapped HA/NA?
        SEGMENT_expected_score=`cat subtype_scores_each_segment.txt | grep \${REF_GENOME_ID} | cut -f \${SEGMENT_POS}`

        SEGMENT_best_counts=`cat subtype_counts_each_segment.txt | sort -rnk\${SEGMENT_POS_COUNTS} | head -1 | cut -f\${SEGMENT_POS_COUNTS}`
        FOUND_SUBTYPES_COUNTS+=(\${SEGMENT_best_counts})

        if [ \${SEGMENT_best} == \${REF_GENOME_ID} ]; then
            # I am not doing anything, I found what I was looking for
            FOUND_SUBTYPES_REASSORTMENT+=(\${REF_GENOME_ID})
            FOUND_SUBTYPES_SCORE_RATIO+=(1)
            FOUND_SUBTYPES_SEQ_SIMILARITY+=(1)

            # Creating a "hybrid" for the genome and BED with primers. We do this N times because
            # if any other segment is not a reassortment, the FASTA files will "grow" with copies
            # of the same sequence.
            cat /home/data/infl/genomes/\${REF_GENOME_ID}/\${REF_GENOME_ID}.fasta | \
                awk -v ID=\${REF_GENOME_ID} '{if (substr(\$0, 1, 1)==">") {filename=("regular_"ID"_"substr(\$0,2) ".fasta");  print \$0 >> filename } else {print toupper(\$0)  >> filename}}'
            cat regular_\${REF_GENOME_ID}_chr?_\${segment}.fasta >> hybrid_genome.fasta
            cat /home/data/infl/primers/\${REF_GENOME_ID}/\${REF_GENOME_ID}_primers.bed | \
                grep \${segment} >> hybrid_primers.bed
            rm regular_\${REF_GENOME_ID}_chr*
        else
            FOUND_SUBTYPES_REASSORTMENT+=(\${SEGMENT_best})
            RATIO=`echo "\${SEGMENT_expected_score} / \${SEGMENT_best_score}" | bc -l`
            FOUND_SUBTYPES_SCORE_RATIO+=(\$RATIO)
            # I calculate the similarity between the segment sequence from the expected subtype
            # and the found one
            cat  /home/data/infl/genomes/\${SEGMENT_best}/\${SEGMENT_best}.fasta | \
                awk -v ID=\${SEGMENT_best} '{if (substr(\$0, 1, 1)==">") {filename=(ID"_"substr(\$0,2) ".fasta");  print \$0"_"ID >> filename } else {print toupper(\$0) >> filename}}'

            # And also a "regular" FASTA for the hybrid
            cat /home/data/infl/genomes/\${SEGMENT_best}/\${SEGMENT_best}.fasta | \
                awk -v ID=\${SEGMENT_best} '{if (substr(\$0, 1, 1)==">") {filename=("regular_"ID"_"substr(\$0,2) ".fasta");  print \$0 >> filename } else {print toupper(\$0) >> filename}}'

            cat \${REF_GENOME_ID}_chr?_\${segment}.fasta \${SEGMENT_best}_chr?_\${segment}.fasta >> tmp.fa
            SEGMENT_alignment_score=`run_nw.py tmp.fa | bc -l`
            FOUND_SUBTYPES_SEQ_SIMILARITY+=(\${SEGMENT_alignment_score=})
            # Remove temporary files for the segment reassortment test
            if awk "BEGIN {exit !(\${SEGMENT_alignment_score} < 0.9 && \${RATIO} < 0.6 && \${SEGMENT_best_counts} >= (${params.min_cov} * 1.5))}"
            then
                echo "Reassortment detected for the segment \${segment}: \${REF_GENOME_ID} -> \${SEGMENT_best}\tAverage coverage is \${SEGMENT_best_counts}"
                cat regular_\${SEGMENT_best}_chr?_\${segment}.fasta >> hybrid_genome.fasta
                cat /home/data/infl/primers/\${SEGMENT_best}/\${SEGMENT_best}_primers.bed | \
                    grep \${segment} >>  hybrid_primers.bed
            elif
                awk "BEGIN {exit !(\${SEGMENT_alignment_score} < 0.9 && \${RATIO} < 0.6 && \${SEGMENT_best_counts} >= (${params.min_cov} * 1.5))}"
            then
                echo "Possible reassortment detected for the segment \${segment}: \${REF_GENOME_ID} -> \${SEGMENT_best}\tAverage coverage is \${SEGMENT_best_counts}"
                cat regular_\${SEGMENT_best}_chr?_\${segment}.fasta >> hybrid_genome.fasta
                cat /home/data/infl/primers/\${SEGMENT_best}/\${SEGMENT_best}_primers.bed | \
                    grep \${segment} >> hybrid_primers.bed
            else
                # The regular FASTA goes into the hybrid.
                cat /home/data/infl/genomes/\${REF_GENOME_ID}/\${REF_GENOME_ID}.fasta | \
                    awk -v ID=\${REF_GENOME_ID} '{if (substr(\$0, 1, 1)==">") {filename=("regular_"ID"_"substr(\$0,2) ".fasta");  print \$0 >> filename } else {print toupper(\$0)  >> filename}}'
                cat regular_\${REF_GENOME_ID}_chr?_\${segment}.fasta >> hybrid_genome.fasta
                cat /home/data/infl/primers/\${REF_GENOME_ID}/\${REF_GENOME_ID}_primers.bed | \
                    grep \${segment} >>  hybrid_primers.bed
                rm regular_\${REF_GENOME_ID}_*
            fi
            rm tmp.fa
            rm *\${SEGMENT_best}_chr*.fasta
        fi
    done

    rm \${REF_GENOME_ID}_chr*.fasta
    echo \${SEGMENTS_REASSORTMENT[@]} >> intermediate.txt
    echo \${FOUND_SUBTYPES_REASSORTMENT[@]} >> intermediate.txt
    echo \${FOUND_SUBTYPES_SCORE_RATIO[@]} >> intermediate.txt
    echo \${FOUND_SUBTYPES_SEQ_SIMILARITY[@]} >> intermediate.txt
    echo \${FOUND_SUBTYPES_COUNTS[@]} >> intermediate.txt

    # That will mess up viewing in IGV...
    REFERENCE_GENOME_FASTA="hybrid_genome.fasta"
    bwa index \${REFERENCE_GENOME_FASTA}
    PRIMERS="hybrid_primers.bed"
    """
}
