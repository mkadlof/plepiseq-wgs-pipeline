// Programs parameters can be modified by a shell wrapper
params.memory = 2024
params.quality_initial = 5
params.min_cov = 20
params.variant = "UNK"

// Specifying location of reads and primers and necessary databases, MUST be selected by a user
params.threads = 12

// output dir, by default "results" directory
params.results_dir = "./results"

include { detect_subtype } from "${params.modules}/infl/detect_subtype.nf"
include { reassortment } from "${params.modules}/infl/reassortment.nf"
include { fastqc as fastqc_1 } from "${params.modules}/common/fastqc.nf"

workflow{
    // Channels
    reads = Channel.fromFilePairs(params.reads)

    // Processes
    fastqc_1(reads, "initialfastq")
    detect_subtype(reads)
    reassortment(detect_subtype.out)
}