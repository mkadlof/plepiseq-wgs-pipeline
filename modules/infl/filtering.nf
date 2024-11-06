process filtering {
    tag "filtering:${sampleId}"
    container  = params.main_image
    input:
    tuple val(sampleId), path(bam), path(bai), path(ref_genome_with_index), val(QC_status), path(primers), path(pairs)
    output:
    tuple val(sampleId), path('to_clip_sorted.bam'), path('to_clip_sorted.bam.bai'), path(primers), path(pairs), val(QC_status), emit: one_amplicon_primers_and_QC
    
    script:
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }
    """
    # Custom filtering for influenza
    # The script takes in sequence:
    # - BAM for filtering and downsampling
    # - Primer scheme
    # - Target coverage per segment (integer)
    # - Minimum read mapping quality (integer)
    # - Reference genome sequences in FASTA (all segments)

    if [ ${QC_status} == "nie" ]; then
      touch to_clip_sorted.bam
      touch to_clip_sorted.bam.bai
    else
      simple_filter_illumina_INFL.py ${bam} ${primers} ${params.max_depth} ${params.min_mapq} ${params.length} ${ref_genome_with_index[final_index]}
    fi
    """
}


process filtering_nanopore {
    tag "filtering:${sampleId}"
    container  = params.main_image
    input:
    tuple val(sampleId), path(bam), path(bai), path(ref_genome_with_index), val(QC_status), path(primers), path(pairs)
    output:
    tuple val(sampleId), path('to_clip_sorted.bam'), path('to_clip_sorted.bam.bai'), path(primers), path(pairs), val(QC_status), emit: one_amplicon_primers_and_QC

    script:
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }
    """


simple_filter_nanopore_final_with_windowstep.py
