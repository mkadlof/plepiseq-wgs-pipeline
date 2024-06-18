// Programs parameters can be modified by a shell wrapper
params.threads = 12
params.variant = "UNK"

include { infl_ref_genome_map } from "${params.modules}/infl/infl_ref_genome_map.nf"


workflow{
    // Channels
    reads = Channel.fromFilePairs(params.reads)

    // Processes
    infl_ref_genome_map(reads)
}