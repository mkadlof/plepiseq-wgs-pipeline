# Module: trimmomatic

## Description  
This module processes raw sequencing reads by trimming adapters and low-quality bases. It uses the Trimmomatic tool to clean the input reads, producing paired and unpaired output files. If the quality control (QC) status is "no," the module creates dummy output files.

## Input  
* tuple val(sampleId), path(reads), val(QC_status)

## Output  
* tuple val(sampleId), path('*_paired.fastq.gz'), val(QC_status), emit: proper_reads_and_qc  
* tuple val(sampleId), path('*_paired.fastq.gz'), emit: proper_reads  
* tuple val(sampleId), path('*_paired.fastq.gz'), path('*_unpaired.fastq.gz'), val(QC_status), emit: all_reads

This module performs quality trimming on sequencing reads based on the QC status and produces cleaned paired and unpaired reads. It emits three types of output: proper reads with QC status, proper reads without QC status, and all reads (paired and unpaired).
