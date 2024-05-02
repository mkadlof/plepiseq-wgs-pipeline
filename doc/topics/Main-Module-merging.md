# Main: Module merging

```Bash
    samtools merge -o clean_sort_dedup_trimmed_sort_tmp.bam ${filtering_bam} ${ivar_bam}
    samtools sort -@ ${params.threads} -o clean_sort_dedup_trimmed_sort.bam clean_sort_dedup_trimmed_sort_tmp.bam
    samtools index clean_sort_dedup_trimmed_sort.bam
```
At this stage, we use Samtools to merge reads belonging to one amplicon with ivar-masked primers from files containing reads from amplicon fusions where one amplicon had low usage. These files also had masked primers.

