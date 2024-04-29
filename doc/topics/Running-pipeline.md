# Running pipeline

## Pipeline parameters

There are three types of flags that we use explicit, implicit and nextflow flags. Explicit MUST be provided during
starting pipeline - they are paths for input files. Explicit parameters are mostly numeric values
for various modules. They have set reasonable default values and usually there is no need to modify them.
NextFlow flags aplay to the way how NextFlow is executed rather than to pipeline itself.

By convention pipeline params starts with two dash `--param_name`, while NextFlow flags starts with single dash `-param-name`.  

### Explicit pipeline parameters

```Bash
    ./run_pipeline.sh \
    --ref_genome 'path/to/reference/genome.fasta' \
    --reads 'path/to/reads/sample_id_{1,2}.fastq.gz' \
    --primers 'path/to/primers.bed' \
    --pairs 'path/to/pairs.tsv' \
    --adapters 'path/to/adapters.fa' \
    --pangolin_db_absolute_path_on_host '/home/user/path/to/pangolin_db' \
    --nextclade_db_absolute_path_on_host '/home/user/path/to/nextclade_db' \
    --kraken2_db_absolute_path_on_host '/home/user/path/to/kraken2_db' \
    --freyja_db_absolute_path_on_host '/home/user/path/to/freyja'
```

`ref_genome` - path to reference genome fasta file

`reads` - path to reads in fastqc format. Must be gzipped and be in form: `sample_id_{1,2}.fastq.gz`. Name must be resorvable by shell into two different files. One for forward reads, and second fo reverse reads.

`primers` - primers in bed file format. Example is below.

```
MN908947.3	2826	2850	nCoV-2019_10_LEFT	1	+	TGAGAAGTGCTCTGCCTATACAGT
MN908947.3	3183	3210	nCoV-2019_10_RIGHT	1	-	TCATCTAACCAATCTTCTTCTTGCTCT
(...)
```

Common primers sets are included in `data/generic/primers` directory and include following:

`SARS1_partmerge_exp`, `SARS2_partmerge_exp`, `V1`, `V2`, `V3`, `V4`, `V4.1`, `V1200`, `V1201`

`pairs` - definition of primers identifiers in two column tab separated file. This file is included in corresponding every primers set in `data/generic/primers`. Structure of primer identifier is meaningful. 
Must match regexp `nCoV-2019_[1,2]_(LEFT,RIGHT)`. Example:

```
nCoV-2019_1_LEFT	nCoV-2019_1_RIGHT
nCoV-2019_2_LEFT	nCoV-2019_2_RIGHT
(...)
```

`adapters` - path to fasta file with adapters. Common adapters are included in `data/generic/primers`.
Example:
```
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

### Implicit pipeline parameters

All implicit parameters are listed in main pipeline file `nf_pipeline.nf` with their reasonable defaults.

```Groovy
params.threads = 5
params.memory = 2024
params.quality_initial = 5
params.length = 90
params.max_number_for_SV = 200000
params.max_depth = 6000
params.min_cov = 20
params.mask = 20
params.quality_snp = 15
params.pval = 0.05
params.lower_ambig = 0.45
params.upper_ambig = 0.55
params.ref_genome_id = "MN908947.3"
```

To run pipeline with modified parameter simply add appropriate flag:

```Bash
./run_nf_pipeline.sh --threads 10
```

`--threads`, positive integer, number of threads used by couple of different modules. In single sample mode should be set to 1/3 of available CPUS. In multisample mode should be adjusted empirically. Recommended value for decent server: 5.

`--memory`, positive integer in MiB, similar to above. Default: 2024

`--quality_initial`, positive integer in PHRED scale, per base quality threshold used in various filtering and reporting modules. Default: 5.

`--length`, positive integer - number of base pairs, minimum length of a read. Default: 90

`--max_number_for_SV` - positive integer, maximum number of reads in bam file for manta module, down sampled by Picard, Default: 200000

`--max_depth` - positive integer, number of base pairs, threshold for short indel callers, reads above this value will be discarded. This is used for speedup indel calling. Default: 6000

`--min_cov` - positive integer, number of base pairs, threshold below which mutation will not be called. Default: 20

`--mask` - positive integer, number of base pairs, below this coverage value, genome will be masked with N. Should be the same as `min_cov`. Default: 20

`--quality_snp` - positive integer, PHRED scale, minimum quality of a base for INDEL calling, Default: 15

`--pval` - float from range [0; 1], minimal probability for INDEL calling, Default: 0.05

`--lower_ambig` and `--upper_ambig` - float from range [0, 1], if fraction of reads introducing alternative allel, fall within this range the position will be classified as *ambiguity*. `upper_ambig` must be greater than `lower_ambig`. Default: [0.45; 0,55]

`--ref_genome_id` - string, identifier of reference genome. Do not change unless you know what you are doing. Default: `MN908947.3`

### Nextflow parameters

This parameters comes with Nextflow and should not be modified without solid reason.

`-config` path to nextflow config file. Default file `nextflow.config` is provided with repo. 
`-with-report` path to report from pipeline execution. May be safely disabled. Default: `report.html`
`-with-dag` path tu file with pipeline graph. May be safely disabled Default: `flowchart-raw.png`
`-with-docker` Docker image used for execution processes. Strictly required. Default: `nf_illumina_sars-3.0-main:latest`
`-resume` Control if restarted pipeline should use cached results or not. Irrelevant in production environment, since every sample will be always run exactly once. In development or during debug may significantly speed up things.


