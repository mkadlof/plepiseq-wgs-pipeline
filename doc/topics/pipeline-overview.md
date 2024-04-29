# Pipeline overview

![flowchart](flowchart.png "Overview of the pipeline.")

NextFlow pipeline for SARS-CoV-2 Illumina data consist of 30 NextFlow modules.
The steps are depicted on image in [](pipeline-overview.md) (`IndexGenome` module which gives input to many other modules was hidden for clarity). Modules are grouped into logical sections, by their function.

### Main section
Goal of **Main** section is to perform mapping with [bwa](%bwa.url%) aligner, filtering, small indel calling by three callers ([varScan](%varscan.url%), [freeBayes](%freebayess.url%) and [lofreq](%lofreq.url%)), identify low quality / coverage regions and finally obtain consensus sequence.

### Quality control 
This section executes [FastQC](%fastqc.url%). Test is run on raw dataset and after trimming adapters with [Trimmomatic](%trimmomatic.url%).

### Contamination detection
This section test if sample is nt contaminated with other species than SarsCov-2.
[Kraken2](%kraken2.url%) is used.

### Dehumanization
This section is used to create a subset of reads in FASTQ format without any reads from organisms other than SARS-CoV-2. Own python scripts are used.

### Coinfections
This section is to detect infections with more than one variant of Sars-Cov-2 virus.
Coinfections are detected with both own python scripts and [Freya](%freyja.url%).

### Functional analysis
This module returns the effects of nucleotide mutations onto protein sequence. We employ [SnpEff](%snpeff.url%).

### SVs Calling
This module is capable to detect large genome rearrangements (*Structural variants*). This is done by down sampling reads with [Picard](%picard.url%) and running [Manta](%manta.url%).

### Statistics
Some simple stats ara calculated with both - own python script and tools from Picard tools.

### Lineages
After SVs calling [pango](%pangolin.url%) line and [nextclade](%nextclade.url%) lineages are identified.

### Building Spike protein model
As a last step Spike protein model is built with [Modeller](%modeller.url%)
