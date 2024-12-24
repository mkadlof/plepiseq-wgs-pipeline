# Main: Module bwa

```Bash
    bwa index ${reference_fasta}
    bwa mem -t ${params.threads} -T 30 ${reference_fasta} ${reads[0]} ${reads[1]} | \
        samtools view -@ ${params.threads} -Sb -f 3 -F 2048 - | \
        samtools sort -@ ${params.threads} -o mapped_reads.bam -
    samtools index mapped_reads.bam
```

We map using the BWA program. For consistency of results, we no longer utilize the option to use a different aligner (such as Bowtie or Minimap2). After mapping, the reads are filtered. We keep only those read pairs in which both reads have been mapped to the reference genome and are proper pairs (-f 3), while simultaneously removing additional/alternative mappings for the reads (-F 2048). The reads are then sorted and indexed.