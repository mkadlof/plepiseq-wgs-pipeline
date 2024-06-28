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
params.window_size = 50 // Window size in which we equalize the coverage
params.min_mapq = 30


// Specifying location of reads and primers and necessary databases, MUST be selected by a user
params.threads = 5
params.reads = "" // Must be provided by user
params.primers_id = "" // We replace the name of the available directory in the wrapper with V4, V4.1, etc. The user provides this.
params.pangolin_db_absolute_path_on_host = ""
params.nextclade_db_absolute_path_on_host = ""
params.kraken2_db_absolute_path_on_host = ""
params.freyja_db_absolute_path_on_host = ""

// adapters, can be set  by a user but there is a default
params.adapters_id="TruSeq3-PE-2" // To podaje user ale jest default

// output dir, by default "results" directory
params.results_dir = "./results"

// Directory with modules, maybe move the .nf file to the main script to remove one settable parameter ? MUST be indicated by user
params.modules = "/absolute/path/to/directory/with/modules"

adapters="/home/data/common/adapters/${params.adapters_id}.fa"

params.ref_genome_id="MN908947.3"

include { kraken2 } from "${params.modules}/sarscov2/kraken2.nf"
include { bwa } from "${params.modules}/common/bwa.nf"
include { dehumanization } from "${params.modules}/sarscov2/dehumanization.nf"
include { fastqc as fastqc_1 } from "${params.modules}/common/fastqc.nf"
include { fastqc as fastqc_2 } from "${params.modules}/common/fastqc.nf"
include { trimmomatic } from "${params.modules}/common/trimmomatic.nf"
include { filtering } from "${params.modules}/sarscov2/filtering.nf"
include { masking } from "${params.modules}/common/masking.nf"
include { merging } from "${params.modules}/sarscov2/merging.nf"
include { picard } from "${params.modules}/common/picard.nf"
include { manta } from "${params.modules}/sarscov2/manta.nf"
include { indelQual } from "${params.modules}/common/indelQual.nf"
include { wgsMetrics } from "${params.modules}/sarscov2/wgsMetrics.nf"
include { lowCov } from "${params.modules}/sarscov2/lowCov.nf"
include { varScan } from "${params.modules}/sarscov2/varscan.nf"
include { freeBayes } from "${params.modules}/sarscov2/freeBayes.nf"
include { lofreq } from "${params.modules}/sarscov2/lofreq.nf"
include { consensus } from "${params.modules}/sarscov2/consensus.nf"
include { vcf_for_fasta } from "${params.modules}/sarscov2/vcf_for_fasta.nf"
include { consensusMasking } from "${params.modules}/sarscov2/consensusMasking.nf"
include { snpEff } from "${params.modules}/sarscov2/snpEff.nf"
include { simpleStats } from "${params.modules}/sarscov2/simpleStats.nf"
include { nextclade } from "${params.modules}/sarscov2/nextclade.nf"
include { pangolin } from "${params.modules}/sarscov2/pangolin.nf"
include { modeller } from "${params.modules}/sarscov2/modeller.nf"

// Coinfection line
include { freyja } from "${params.modules}/sarscov2/freyja.nf"
include { coinfection_ivar } from "${params.modules}/sarscov2/coinfection_ivar.nf"
include { coinfection_varscan } from "${params.modules}/sarscov2/coinfection_varscan.nf"
include { coinfection_analysis } from "${params.modules}/sarscov2/coinfection_analysis.nf"

workflow{
    // Channels
    reads = Channel.fromFilePairs(params.reads)
    ref_genome = Channel.fromFilePairs("${projectDir}/data/sarscov2/genome/*sarscov2.fasta", size: -1).collect().sort()
    ref_genome_with_index = Channel.fromFilePairs("${projectDir}/data/sarscov2/genome/*sarscov2.fasta{,.amb,.ann,.bwt,.fai,.pac,.sa}", size: -1).collect().sort()
    primers_and_pairs = Channel.fromFilePairs("${projectDir}/data/sarscov2/primers/${params.primers_id}/*{nCoV-2019.scheme.bed,pairs.tsv}").collect()
    coinfections = Channel.fromPath("${projectDir}/data/sarscov2/coinfections/").first()
    vcf_template = Channel.fromPath("${projectDir}/data/sarscov2/vcf_template/").first()
    modeller_data = Channel.fromPath("${projectDir}/data/sarscov2/modeller/").first()

    // The following construct is needed to add the sampleId to the channels
    // with constant data (like reference genome or primers). This is done to
    // generalize cases where the reference genome may not be identical for
    // each sample (as is the case in the influenza pipeline, for example).
    ref_genome_with_index = ref_genome_with_index
    .combine(reads).map { _, genome_and_index_files, sampleId, reads ->
        return [sampleId, genome_and_index_files]
    }
    primers_and_pairs = primers_and_pairs
    .combine(reads).map { _, primers_and_pairs, sampleId, reads ->
        return [sampleId, primers_and_pairs]
    }

    // Processes
    fastqc_1(reads, "initialfastq")
    kraken2(reads)
    trimmomatic(reads, adapters)
    j = trimmomatic.out[0].join(ref_genome_with_index)
    bwa(j)
    dehumanization(bwa.out, trimmomatic.out[1])
    fastqc_2(trimmomatic.out[0], "aftertrimmomatic")
    filtering(bwa.out, primers_and_pairs.map{sampleId, files -> [sampleId, files[0]]})
    masking(filtering.out[0], primers_and_pairs)
    combined = filtering.out[1].join(masking.out)
    merging(combined)
    picard(bwa.out)
    indelQual(merging.out, ref_genome)
    wgsMetrics(indelQual.out, ref_genome)
    lowCov(indelQual.out, ref_genome)
    varScan(indelQual.out, ref_genome)
    freeBayes(indelQual.out, ref_genome)
    lofreq(indelQual.out, ref_genome_with_index)
    c1 = varScan.out.join(freeBayes.out).join(lofreq.out)
    consensus(c1)
    vcf_for_fasta(consensus.out, ref_genome, vcf_template)
    consensusMasking(consensus.out.join(lowCov.out[1]))
    manta(picard.out.join(consensusMasking.out))
    nextclade(manta.out)
    modeller(nextclade.out[1], modeller_data)
    pangolin(manta.out)
    snpEff(vcf_for_fasta.out.join(indelQual.out), ref_genome)
    simpleStats(manta.out.join(wgsMetrics.out), primers_and_pairs.map{sampleId, files -> [sampleId, files[0]]})

    // Coinfection line
    coinfection_ivar(bwa.out, ref_genome, primers_and_pairs.map{sampleId, files -> [sampleId, files[0]]})
    freyja(coinfection_ivar.out[0], ref_genome)
    coinfection_varscan(coinfection_ivar.out[1])
    coinfection_analysis(coinfection_varscan.out, coinfections)
}
