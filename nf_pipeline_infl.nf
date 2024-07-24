// Programs parameters can be modified by a shell wrapper
params.memory = 2024
params.quality_initial = 5
params.length = 90
params.max_number_for_SV = 200000
params.max_depth = 600
params.min_cov = 20
params.mask = 20
params.quality_snp = 15
params.pval = 0.05
params.lower_ambig = 0.45
params.upper_ambig = 0.55
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
nextalign_db="/home/data/infl/nextalign"

include { detect_subtype } from "${params.modules}/infl/detect_subtype.nf"
include { reassortment } from "${params.modules}/infl/reassortment.nf"
include { bwa } from "${params.modules}/common/bwa.nf"
include { fastqc as fastqc_1 } from "${params.modules}/common/fastqc.nf"
include { fastqc as fastqc_2 } from "${params.modules}/common/fastqc.nf"
include { trimmomatic } from "${params.modules}/common/trimmomatic.nf"
include { filtering } from "${params.modules}/infl/filtering.nf"
include { masking } from "${params.modules}/common/masking.nf"
include { picard } from "${params.modules}/common/picard.nf"
include { sort_and_index } from "${params.modules}/infl/sort_and_index.nf"
include { indelQual } from "${params.modules}/common/indelQual.nf"
include { wgsMetrics } from "${params.modules}/common/wgsMetrics.nf"
include { lowCov } from "${params.modules}/common/lowCov.nf"
include { varScan } from "${params.modules}/common/varscan.nf"
include { freeBayes } from "${params.modules}/common/freeBayes.nf"
include { lofreq } from "${params.modules}/common/lofreq.nf"
include { consensus } from "${params.modules}/common/consensus.nf"
include { nextclade } from "${params.modules}/infl/nextclade.nf"
include { nextalign } from "${params.modules}/infl/nextalign.nf"

workflow{
    // Channels
    reads = Channel.fromFilePairs(params.reads)
    genomes = Channel.fromPath("${projectDir}/data/infl/genomes/")
    primers = Channel.fromPath("${projectDir}/data/infl/primers/")
    pairs = Channel.fromPath("${projectDir}/data/infl/primers/pairs.tsv").first()

    // Processes
    fastqc_1(reads, "initialfastq")
    trimmomatic(reads, adapters)
    fastqc_2(trimmomatic.out[0], "aftertrimmomatic")
    detect_subtype(reads, genomes)
    reassortment(detect_subtype.out[0], genomes, primers)
    // For convenience we name the output with the hybrid genome as ref_genome
    ref_genome = reassortment.out[0].map{ sampleId, files -> [sampleId, files[0]]}
    bwa(trimmomatic.out[0].join(reassortment.out[0]))
    c1 = bwa.out.join(reassortment.out[0]).join(reassortment.out[1])
    filtering(c1)

    primers_and_pairs = reassortment.out[1].merge(pairs).map {sampleId, primers, pairs ->
        return [sampleId, tuple(primers, pairs)]}

    masking(filtering.out, primers_and_pairs)
    picard(bwa.out)
    sort_and_index(masking.out)
    indelQual(sort_and_index.out, ref_genome)
    varScan(indelQual.out, ref_genome)
    freeBayes(indelQual.out, ref_genome)
    lofreq(indelQual.out, ref_genome)
    wgsMetrics(indelQual.out, ref_genome)
    lowCov(indelQual.out, ref_genome)
    c2 = lowCov.out[1].join(varScan.out).join(freeBayes.out).join(lofreq.out)
    consensus(c2)
    nextclade(detect_subtype.out[1], consensus.out[1])
    nextalign(detect_subtype.out[1], consensus.out[1], nextalign_db)
}