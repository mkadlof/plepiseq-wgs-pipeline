// Number of threads (used in several different modules)
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
params.quality_SNP = 15
params.ref_genome_id = "MN908947.3"


include { indexGenome } from './modules/indexGenome.nf'
include { bwa } from './modules/bwa.nf'
include { fastqc as fastqc_1 } from './modules/fastqc.nf'
include { fastqc as fastqc_2 } from './modules/fastqc.nf'
include { trimmomatic } from './modules/trimmomatic.nf'
include { filtering } from './modules/filtering.nf'
include { masking } from './modules/masking.nf'
include { merging } from './modules/merging.nf'
include { picard } from './modules/picard.nf'
include { viterbi } from './modules/viterbi.nf'
include { wgsMetrics } from './modules/wgsMetrics.nf'
include { lowCov } from './modules/lowCov.nf'
include { varScan } from './modules/varscan.nf'
include { freeBayes } from './modules/freeBayes.nf'
include { lofreq } from './modules/lofreq.nf'
include { consensus } from './modules/consensus.nf'
include { consensusMasking } from './modules/consensusMasking.nf'
include { variantIdentification } from './modules/variantIdentification.nf'
include { functionalAnalysis } from './modules/functionalAnalysis.nf'
include { consensusAnalysis } from './modules/consensusAnalysis.nf'


workflow{
    // Channels
    ref_genome = Channel.value(params.ref_genome as Path)
    reads = Channel.fromFilePairs(params.reads)
    primers = Channel.value(params.primers as Path)
    pairs = Channel.value(params.pairs as Path)
    adapters = Channel.value(params.adapters as Path)

    // Processes
    indexGenome(ref_genome)
    fastqc_1(reads, "initialfastq")
    trimmomatic(reads, adapters)
    bwa(trimmomatic.out[0], indexGenome.out)
    fastqc_2(trimmomatic.out[0], "aftertrimmomatic")
    filtering(bwa.out, primers)
    masking(filtering.out[0], primers, pairs)
    combined = filtering.out[1].join(masking.out)
    merging(combined, primers, pairs)
    picard(bwa.out)
    viterbi(merging.out, indexGenome.out)
    wgsMetrics(viterbi.out, indexGenome.out)
    lowCov(viterbi.out, indexGenome.out)
    varScan(viterbi.out.join(lowCov.out[1]), indexGenome.out)
    freeBayes(viterbi.out.join(lowCov.out[1]), indexGenome.out)
    lofreq(viterbi.out.join(lowCov.out[1]), indexGenome.out)
    consensus(varScan.out[1].join(freeBayes.out[1]).join(lofreq.out[1]))
    consensusMasking(consensus.out.join(lowCov.out[1]))
    variantIdentification(varScan.out[1].join(freeBayes.out[1]).join(lofreq.out[1]).join(consensusMasking.out))
    functionalAnalysis(varScan.out[0].join(freeBayes.out[0]).join(lofreq.out[0]))
    consensusAnalysis(varScan.out[0].join(freeBayes.out[0]).join(lofreq.out[0]))
}