process coinfection_genome_masking_illumina {
    tag "coinfection_genome_masking:${sampleId}"
    container  = params.main_image

    input:
    tuple val(sampleId), path(mapped_reads), path(mapped_reads_bai), path(ref_genome_with_index), path(primers), val(QC_status)

    output:
    tuple val(sampleId), path('for_contamination_sorted.bam'), path('for_contamination_sorted.bam.bai'), path(ref_genome_with_index), val(QC_status), emit: to_freyja
    tuple val(sampleId), path('for_contamination.mpileup'), val(QC_status), emit: to_custom_analysis

    script:
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }
    """
    if [ ${QC_status} == "nie" ]; then
      touch for_contamination_sorted.bam for_contamination_sorted.bam.bai
      touch for_contamination.mpileup
    else
      ivar trim -i ${mapped_reads} \
                -b ${primers} \
                -m ${params.length} \
                -q ${params.quality_initial} \
                -e \
                -p for_contamination

      samtools sort -@ ${params.threads} -o for_contamination_sorted.bam for_contamination.bam
      samtools index for_contamination_sorted.bam
      samtools mpileup --max-depth 10000 \
                   --fasta-ref  ${ref_genome_with_index[final_index]} \
                   --min-BQ ${params.quality_snp} \
                   for_contamination_sorted.bam >> for_contamination.mpileup
    fi
    """
}

process coinfection_genome_masking_nanopore {
    tag "coinfection_genome_masking:${sampleId}"
    container  = params.main_image

    input:
    tuple val(sampleId), path(bam), path(bai), path(ref_genome_with_index), path(primers), val(QC_status)

    output:
    tuple val(sampleId), path('for_contamination_sorted.bam'), path('for_contamination_sorted.bam.bai'), path(ref_genome_with_index), val(QC_status), emit: to_freyja
    tuple val(sampleId), path('for_contamination.mpileup'), val(QC_status), emit: to_custom_analysis

    script:
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }
    """
    if [ ${QC_status} == "nie" ]; then
      touch for_contamination_sorted.bam for_contamination_sorted.bam.bai
      touch for_contamination.mpileup
    else
      
    samtools ampliconclip --filter-len 1 --both-ends -b ${primers} --tolerance ${params.bed_offset} -o for_contaminations_presorted.bam -O bam ${bam}
    samtools sort  -@ ${params.threads} -o for_contamination_sorted.bam for_contaminations_presorted.bam
    samtools index for_contamination_sorted.bam
    samtools mpileup  -B -f ${ref_genome_with_index[final_index]} -Q 1 for_contamination_sorted.bam >> for_contamination.mpileup

    fi
    """
}

