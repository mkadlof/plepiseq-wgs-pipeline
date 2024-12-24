# Module: copy_genome_and_primers

## Description  
This is an auxiliary module aimed at enabling the use of the same modules for SARS as for influenza, for which a hybrid reference genome is being built. This module manages the copying of genome and primer data into the Nextflow working directory. It handles both Illumina and Nanopore data by copying the relevant files depending on the QC status. If QC fails, dummy files are generated.

## Input  
* tuple val(sampleId), path(reads), val(QC_STATUS)

## Output  
* tuple val(sampleId), path("sars*"), path("primers.bed"), path("pairs.tsv"), env(REF_GENOME_ID), env(QC_exit), emit: all  
* tuple val(sampleId), path("genome.fasta"), path("primers.bed"), env(QC_exit), emit: all_nanopore  
* tuple val(sampleId), path("sars*"), path("primers.bed"), env(QC_exit), emit: to_bwa  
* tuple val(sampleId), path("primers.bed"), path("pairs.tsv"), emit: primers_and_pairs  
* tuple val(sampleId), path("primers.bed"), emit: primers  
* tuple val(sampleId), path("genome.fasta"), emit: only_genome

This module copies genome and primer files to the work directory, preparing the data for downstream analysis. It checks the QC status and either copies real data or generates dummy files when QC fails. It also sets environment variables related to the reference genome ID and QC status.
