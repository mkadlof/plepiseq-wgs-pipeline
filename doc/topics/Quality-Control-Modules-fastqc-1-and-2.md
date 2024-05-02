# Quality Control: Modules fastqc 1 and 2

```Bash
    fastqc --format fastq \
           --threads ${params.threads} \
           --memory ${params.memory} \
           --extract \
           --delete \
           --outdir . \
           ${reads[0]} ${reads[1]}
    r1=\$(basename ${reads[0]} .fastq.gz)
    r2=\$(basename ${reads[1]} .fastq.gz)
    parse_fastqc_output.py \${r1}_fastqc/fastqc_data.txt ${params.quality_initial} >> ${prefix}_forward_fastqc.txt
    parse_fastqc_output.py \${r2}_fastqc/fastqc_data.txt ${params.quality_initial} >> ${prefix}_reverse_fastqc.txt
  
```

FastQC modules are used for general sequencing quality assessment and detection of common sample issues. The module is run twice, first for the raw reads provided by the user, and then after the action of the trimmomatic program. The result is parsed by a custom Python script, and the results are returned to the users, allowing them to independently assess the sample.
