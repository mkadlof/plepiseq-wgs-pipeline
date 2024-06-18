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
params.mapping_quality = 30


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
params.modules = "/absolute/path/to/directory/with/modules/sarscov2"

// Old parameter not settable anymore by a shell wrapper
params.primers="/home/data/sarscov2/primers/${params.primers_id}/nCoV-2019.scheme.bed"
params.pairs="/home/data/sarscov2/primers/${params.primers_id}/pairs.tsv"
params.adapters="/home/data/sarscov2/adapters/${params.adapters_id}.fa"

params.ref_genome_id="MN908947.3"

include { kraken2 } from "${params.modules}/kraken2.nf"
include { bwa } from "${params.modules}/bwa.nf"
include { dehumanization } from "${params.modules}/dehumanization.nf"
include { fastqc as fastqc_1 } from "${params.modules}/fastqc.nf"
include { fastqc as fastqc_2 } from "${params.modules}/fastqc.nf"
include { trimmomatic } from "${params.modules}/trimmomatic.nf"
include { filtering } from "${params.modules}/filtering.nf"
include { masking } from "${params.modules}/masking.nf"
include { merging } from "${params.modules}/merging.nf"
include { picard } from "${params.modules}/picard.nf"
include { manta } from "${params.modules}/manta.nf"
include { indelQual } from "${params.modules}/indelQual.nf"
include { wgsMetrics } from "${params.modules}/wgsMetrics.nf"
include { lowCov } from "${params.modules}/lowCov.nf"
include { varScan } from "${params.modules}/varscan.nf"
include { freeBayes } from "${params.modules}/freeBayes.nf"
include { lofreq } from "${params.modules}/lofreq.nf"
include { consensus } from "${params.modules}/consensus.nf"
include { vcf_for_fasta } from "${params.modules}/vcf_for_fasta.nf"
include { consensusMasking } from "${params.modules}/consensusMasking.nf"
include { snpEff } from "${params.modules}/snpEff.nf"
include { simpleStats } from "${params.modules}/simpleStats.nf"
include { nextclade } from "${params.modules}/nextclade.nf"
include { pangolin } from "${params.modules}/pangolin.nf"
include { modeller } from "${params.modules}/modeller.nf"

// Coinfection line
include { freyja } from "${params.modules}/freyja.nf"
include { coinfection_ivar } from "${params.modules}/coinfection_ivar.nf"
include { coinfection_varscan } from "${params.modules}/coinfection_varscan.nf"
include { coinfection_analysis } from "${params.modules}/coinfection_analysis.nf"

workflow{
    // Channels
    reads = Channel.fromFilePairs(params.reads)

    // Processes
    fastqc_1(reads, "initialfastq")
    kraken2(reads)
    trimmomatic(reads, params.adapters)
    bwa(trimmomatic.out[0])
    dehumanization(bwa.out, trimmomatic.out[1])
    fastqc_2(trimmomatic.out[0], "aftertrimmomatic")
    filtering(bwa.out, params.primers)
    masking(filtering.out[0], params.primers, params.pairs)
    combined = filtering.out[1].join(masking.out)
    merging(combined)
    picard(bwa.out)
    indelQual(merging.out)
    wgsMetrics(indelQual.out)
    lowCov(indelQual.out)
    varScan(indelQual.out)
    freeBayes(indelQual.out)
    lofreq(indelQual.out)
    c1 = varScan.out.join(freeBayes.out).join(lofreq.out)
    consensus(c1)
    vcf_for_fasta(consensus.out)
    consensusMasking(consensus.out.join(lowCov.out[1]))
    manta(picard.out.join(consensusMasking.out))
    nextclade(manta.out)
    modeller(nextclade.out[1])
    pangolin(manta.out)
    snpEff(vcf_for_fasta.out.join(indelQual.out))
    simpleStats(manta.out.join(wgsMetrics.out))

    // Coinfection line
    coinfection_ivar(bwa.out, params.primers)
    freyja(coinfection_ivar.out[0])
    coinfection_varscan(coinfection_ivar.out[1])
    coinfection_analysis(coinfection_varscan.out)
}
