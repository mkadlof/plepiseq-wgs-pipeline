# Main: Module viterbi

```Bash
    lofreq viterbi --ref ${reference_fasta} \
                   --out clean_sort_dedup_trimmed_sort_viterbi.bam \
                   ${bam}

    lofreq indelqual --ref ${reference_fasta} \
                     --out forvariants.bam \
                     --dindel clean_sort_dedup_trimmed_sort_viterbi.bam
    samtools sort -@ ${params.threads} \
                  -o forvariants.bam \
                  forvariants.bam
    samtools index forvariants.bam
```

This fragment is part of the procedure recommended for identifying reads using the Lofreq program following GATK recommendations. Reads are realigned and a quality score is introduced for indel quality, enabling their identification by the Lofreq program. The above procedure does not affect the results of Varscan or Freebayes.
