// Number of threads (used in several different modules)
params.threads = 5
params.memory = 2024
params.quality_initial = 5

include { bwa } from './modules/bwa.nf'
include { fastqc as fastqc_1 } from './modules/fastqc.nf'

workflow{
    // Channels
    reference_genome = Channel.value(params.reference_genome as Path)
    reads = Channel.fromFilePairs(params.reads)
    primers = Channel.value(params.primers as Path)

    // Processes
    fastqc_1(reads)
    bwa(reads, reference_genome)
}