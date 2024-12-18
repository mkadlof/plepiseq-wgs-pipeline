# Pipeline overview

The pipeline was written using the Nextflow Framework. It has a modular structure in the sense that individual processes (tools) are located in separate files called modules, and the entire pipeline is invoked from the main file named nf_pipeline_viral.nf.

The main file actually contains 6 pipelines for all combinations of three viral organisms (SARS-CoV-2, Influenza, and RSV), and two technologies (Illumina and Nanopore). The pipelines can share selected modules or groups of modules, but in areas requiring specific analyses, they use dedicated modules designed only for those purposes. Which pipeline is invoked is determined by the input parameters.

In each pipeline, the following stages of analysis can be distinguished:

1. Quality control of reads (including contamination detection)
2. Mapping of reads to the reference genome
3. Filtering and downsampling of reads
4. Detection of structural variants
5. Functional analysis
6. Descriptive statistics
7. Phylogenetic analysis
8. 3D Modeling of selected proteins
9. Results aggregation

All pipelines are embedded in Docker images. Containers are run by natively installed Nextflow Framework. 

## Pipelines as DAGs

Each Nextflow module consists of four basic sections:

 - Module settings (meta-information, as well as the method of execution and result return)
 - Definition of input channels
 - Definition of output channels
 - Script executing the task (process)

Taka organizacja pozwala na łączenie modułów w dowolnych kombinacjach, o ile wyjście jednego modułu jest kompatybilne z wyjściem drugiego modułu. Połączone moduły tworzą acykliczny graf skierowany (tzw. DAG). Grafy te wygodnie jest przedstawiać w formie graficznej, dzięki czemu można z łatwo prześledzić proces analizy danych. Poniższe 6 stron zawiera wszystkie 6 pipelineów w formie DAG'ów. Na przedstawionych grafikach owale reprezentują procesy, zaś linie je łączące kanały przepływu danych.

<img alt="Illumina Sars-CoV-2 Flowchart" src="illumina_sars.png" title="Illumina Sars-CoV-2 Flowchart" height="900"/>
<img alt="Illumina Influenza Flowchart" src="illumina_infl.png" title="Illumina Influenza Flowchart" height="900"/>
<img alt="Illumina RSV Flowchart" src="illumina_rsv.png" title="Illumina RSV Flowchart" height="900"/>
<img alt="Nanopore Sars-CoV-2 Flowchart" src="nanopore_sars.png" title="Nanopore Sars-CoV-2 Flowchart" height="900"/>
<img alt="Nanopore Influenza Flowchart" src="nanopore_infl.png" title="Nanopore Influenza Flowchart" height="900"/>
<img alt="Nanopore RSV Flowchart" src="nanopore_rsv.png" title="Nanopore RSV Flowchart" height="900"/>

[//]: # (### Main section)

[//]: # (Goal of **Main** section is to perform mapping with [bwa]&#40;%bwa.url%&#41; aligner, filtering, small indel calling by three callers &#40;[varScan]&#40;%varscan.url%&#41;, [freeBayes]&#40;%freebayess.url%&#41; and [lofreq]&#40;%lofreq.url%&#41;&#41;, identify low quality / coverage regions and finally obtain consensus sequence.)

[//]: # ()
[//]: # (### Quality control )

[//]: # (This section executes [FastQC]&#40;%fastqc.url%&#41;. Test is run on raw dataset and after trimming adapters with [Trimmomatic]&#40;%trimmomatic.url%&#41;.)

[//]: # ()
[//]: # (### Contamination detection)

[//]: # (This section test if sample is nt contaminated with other species than SarsCov-2.)

[//]: # ([Kraken2]&#40;%kraken2.url%&#41; is used.)

[//]: # ()
[//]: # (### Dehumanization)

[//]: # (This section is used to create a subset of reads in FASTQ format without any reads from organisms other than SARS-CoV-2. Own python scripts are used.)

[//]: # ()
[//]: # (### Coinfections)

[//]: # (This section is to detect infections with more than one variant of Sars-Cov-2 virus.)

[//]: # (Coinfections are detected with both own python scripts and [Freya]&#40;%freyja.url%&#41;.)

[//]: # ()
[//]: # (### Functional analysis)

[//]: # (This module returns the effects of nucleotide mutations onto protein sequence. We employ [SnpEff]&#40;%snpeff.url%&#41;.)

[//]: # ()
[//]: # (### SVs Calling)

[//]: # (This module is capable to detect large genome rearrangements &#40;*Structural variants*&#41;. This is done by down sampling reads with [Picard]&#40;%picard.url%&#41; and running [Manta]&#40;%manta.url%&#41;.)

[//]: # ()
[//]: # (### Statistics)

[//]: # (Some simple stats ara calculated with both - own python script and tools from Picard tools.)

[//]: # ()
[//]: # (### Lineages)

[//]: # (After SVs calling [pango]&#40;%pangolin.url%&#41; line and [nextclade]&#40;%nextclade.url%&#41; lineages are identified.)

[//]: # ()
[//]: # (### Building Spike protein model)

[//]: # (As a last step Spike protein model is built with [Modeller]&#40;%modeller.url%&#41;)
