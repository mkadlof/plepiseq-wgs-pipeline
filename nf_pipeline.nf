// Number of threads (used in several different modules)
params.threads = 5

include { bwa } from './modules/bwa.nf'

workflow{
    // Channels
    reference_genome = Channel.value(params.reference_genome as Path)
    reads = Channel.fromFilePairs(params.reads)
    primers = Channel.value(params.primers as Path)
    bwa(reads, reference_genome)
}