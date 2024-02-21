// Number of threads (used in several different modules)
params.threads = 5
params.memory = 2024
params.quality_initial = 5
params.length = 90

include { bwa } from './modules/bwa.nf'
include { fastqc as fastqc_1 } from './modules/fastqc.nf'
include { fastqc as fastqc_2 } from './modules/fastqc.nf'
include { trimmomatic } from './modules/trimmomatic.nf'
include { filtering } from './modules/filtering.nf'
include { masking } from './modules/masking.nf'


workflow{
    // Channels
    reference_genome = Channel.value(params.reference_genome as Path)
    reads = Channel.fromFilePairs(params.reads)
    primers = Channel.value(params.primers as Path)
    pairs = Channel.value(params.pairs as Path)
    adapters = Channel.value(params.adapters as Path)

    // Processes
    fastqc_1(reads, "initialfastq")
    bwa(reads, reference_genome)
    trimmomatic(reads, adapters)
    fastqc_2_input = trimmomatic.out.map { sampleId, forward_paired, forward_unpaired, reverse_paired, reverse_unpaired -> [sampleId, [forward_paired, reverse_paired]] }
    fastqc_2(fastqc_2_input, "aftertrimmomatic")
    filtering(bwa.out, primers)
    masking_input = filtering.out.map { sampleId, bam1, bai1, bam2, bai2, stats -> [sampleId, bam1, bai1] }
    masking(masking_input, primers, pairs)
}