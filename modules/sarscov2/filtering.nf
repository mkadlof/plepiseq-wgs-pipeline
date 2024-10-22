process filtering {
    tag "filtering:${sampleId}"

    input:
    tuple val(sampleId), path(bam), path(bai), val(QC_status), path(primers), path(pairs)

    output:
    tuple val(sampleId), path('first_pass_sorted.bam'), path('first_pass_sorted.bam.bai'), path(primers), path(pairs), val(QC_status), emit: one_amplicon_primers_and_QC
    tuple val(sampleId), path('two_amplicons_sorted.bam'), emit: two_amplicon_only
    // tuple val(sampleId), path('Statystyki.txt')

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch first_pass_sorted.bam
      touch first_pass_sorted.bam.bai
      touch two_amplicons_sorted.bam
      QC_exit="nie"
    else
      bed_offset=0
      length=`echo "${params.length} - 20" | bc -l`
      equal_depth=`echo "${params.max_depth} / 2" | bc -l | awk '{print int(\$0)}'` 
      simple_filter_illumina_one_segment.py ${bam} ${primers} \${bed_offset} \${length} ${params.min_mapq} ${params.window_size} \${equal_depth}
      samtools index first_pass_sorted.bam
    
      if [ -e two_amplicons_sorted.bam ]; then
          cp two_amplicons_sorted.bam tmp.bam
      else
          samtools view -b -o two_amplicons_sorted.bam -L /dev/null first_pass_sorted.bam
      fi
    fi
    """
}
