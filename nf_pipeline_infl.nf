// Programs parameters can be modified by a shell wrapper
params.threads = 12
params.variant = "UNK"
params.min_cov = 20

include { detect_subtype } from "${params.modules}/infl/detect_subtype.nf"
include { reassortment } from "${params.modules}/infl/reassortment.nf"


workflow{
    // Channels
    reads = Channel.fromFilePairs(params.reads)

    // Processes
    detect_subtype(reads)
    reassortment(detect_subtype.out)
}