// Programs parameters can be modified by a shell wrapper
params.memory = 2024
params.quality_initial = 5
params.length = 90
params.max_depth = 600
params.min_cov = 20
params.min_mapq = 30
params.variant = "UNK"

// Specifying location of reads and primers and necessary databases, MUST be selected by a user
params.threads = 12

// adapters, can be set  by a user but there is a default
params.adapters_id="TruSeq3-PE-2" // To podaje user ale jest default

// output dir, by default "results" directory
params.results_dir = "./results"

// Translate id into paths
adapters="/home/data/common/adapters/${params.adapters_id}.fa"


include { detect_subtype } from "${params.modules}/infl/detect_subtype.nf"
include { reassortment } from "${params.modules}/infl/reassortment.nf"
include { bwa } from "${params.modules}/common/bwa.nf"
include { fastqc as fastqc_1 } from "${params.modules}/common/fastqc.nf"
include { fastqc as fastqc_2 } from "${params.modules}/common/fastqc.nf"
include { trimmomatic } from "${params.modules}/common/trimmomatic.nf"
include { filtering } from "${params.modules}/infl/filtering.nf"
include { masking } from "${params.modules}/common/masking.nf"

workflow{
    // Channels
    reads = Channel.fromFilePairs(params.reads)
    genomes = Channel.fromPath("${projectDir}/data/infl/genomes/")
    primers = Channel.fromPath("${projectDir}/data/infl/primers/")

    // Processes
    fastqc_1(reads, "initialfastq")
    trimmomatic(reads, adapters)
    fastqc_2(trimmomatic.out[0], "aftertrimmomatic")
    detect_subtype(reads, genomes)
    reassortment(detect_subtype.out[0], genomes, primers)
    bwa(trimmomatic.out[0], reassortment.out[0])
    filtering(bwa.out, reassortment.out)
    masking(filtering.out, reassortment.out[0], reassortment.out[1])
}