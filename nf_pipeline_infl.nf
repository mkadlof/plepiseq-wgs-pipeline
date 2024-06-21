// Programs parameters can be modified by a shell wrapper
params.memory = 2024
params.quality_initial = 5
params.length = 90
params.min_cov = 20
params.variant = "UNK"

// Specifying location of reads and primers and necessary databases, MUST be selected by a user
params.threads = 12

// adapters, can be set  by a user but there is a default
params.adapters_id="TruSeq3-PE-2" // To podaje user ale jest default

// output dir, by default "results" directory
params.results_dir = "./results"

// Old parameter not settable anymore by a shell wrapper
params.adapters="/home/data/common/adapters/${params.adapters_id}.fa"

include { detect_subtype } from "${params.modules}/infl/detect_subtype.nf"
include { reassortment } from "${params.modules}/infl/reassortment.nf"
include { bwa } from "${params.modules}/common/bwa.nf"
include { fastqc as fastqc_1 } from "${params.modules}/common/fastqc.nf"
include { fastqc as fastqc_2 } from "${params.modules}/common/fastqc.nf"
include { trimmomatic } from "${params.modules}/common/trimmomatic.nf"

workflow{
    // Channels
    reads = Channel.fromFilePairs(params.reads)

    // Processes
    fastqc_1(reads, "initialfastq")
    trimmomatic(reads, params.adapters)
    bwa(trimmomatic.out[0])
    fastqc_2(trimmomatic.out[0], "aftertrimmomatic")
    detect_subtype(reads)
    reassortment(detect_subtype.out)
}