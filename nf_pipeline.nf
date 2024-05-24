// Variables not present in run_nf_pipeline.sh.template file
// See run_nf_pipeline.sh.template menu help for explanations 
params.memory = 2024
params.quality_initial = 5
params.length = 90
params.max_number_for_SV = 200000
params.max_depth = 600 // updated to better fit new filtering approach
params.min_cov = 20
params.mask = 20
params.quality_snp = 15
params.pval = 0.05
params.lower_ambig = 0.45
params.upper_ambig = 0.55
//params.ref_genome_id = "MN908947.3"
params.window_size = 50 // Wielkosc okna w ktorym wyrownujemy pokrycie
params.mapping_quality = 30


// Variables that we strongly recommed to be set via run_nf_pipeline.sh.template file
// already present in that file


// Specifing location of reads/primers/adapters etc.
params.threads = 5
params.ref_genome = ""
params.reads = ""
params.primers = ""
params.pairs = ""
params.adapters = ""


// Subcategory
// paths to EXTERNAL databases with no defaults
params.pangolin_db_absolute_path_on_host = ""
params.nextclade_db_absolute_path_on_host = ""
params.kraken2_db_absolute_path_on_host = ""
params.freyja_db_absolute_path_on_host = ""

params.modules = "/home/michall/git/nf_illumina_sars_ml/modules"
params.results_dir = "./results"



// Old parameter not settable by sh wrapper
params.ref_genome_id=new File(params.ref_genome).readLines().get(0).replaceAll('>','')


include { indexGenome } from "${params.modules}/indexGenome.nf"
include { kraken2 } from "${params.modules}/kraken2.nf"
include { bwa } from "${params.modules}/bwa.nf"
include { dehumanization } from "${params.modules}/dehumanization.nf"
include { fastqc as fastqc_1 } from "${params.modules}/fastqc.nf"
include { fastqc as fastqc_2 } from "${params.modules}/fastqc.nf"
include { trimmomatic } from "${params.modules}/trimmomatic.nf"
include { filtering } from "${params.modules}/filtering_novel.nf"
include { masking } from "${params.modules}/masking.nf"
include { merging } from "${params.modules}/merging.nf"
include { picard } from "${params.modules}/picard.nf"
include { manta } from "${params.modules}/manta.nf"
include { viterbi } from "${params.modules}/viterbi.nf"
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
    ref_genome = Channel.value(params.ref_genome as Path)
    reads = Channel.fromFilePairs(params.reads)
    primers = Channel.value(params.primers as Path)
    pairs = Channel.value(params.pairs as Path)
    adapters = Channel.value(params.adapters as Path)

    // Processes
    indexGenome(ref_genome)
    fastqc_1(reads, "initialfastq")
    kraken2(reads)
    trimmomatic(reads, adapters)
    bwa(trimmomatic.out[0], indexGenome.out)
    dehumanization(bwa.out, trimmomatic.out[1])
    fastqc_2(trimmomatic.out[0], "aftertrimmomatic")
    filtering(bwa.out, primers)
    masking(filtering.out[0], primers, pairs)
    combined = filtering.out[1].join(masking.out)
    merging(combined, primers, pairs)
    picard(bwa.out)
    viterbi(merging.out, indexGenome.out)
    wgsMetrics(viterbi.out, indexGenome.out)
    lowCov(viterbi.out, indexGenome.out)
    varScan(viterbi.out, indexGenome.out)
    freeBayes(viterbi.out, indexGenome.out)
    lofreq(viterbi.out, indexGenome.out)
    c1 = varScan.out.join(freeBayes.out).join(lofreq.out)
    consensus(c1)
    vcf_for_fasta(consensus.out, indexGenome.out)
    consensusMasking(consensus.out.join(lowCov.out[1]))
    manta(picard.out.join(consensusMasking.out), indexGenome.out)
    nextclade(manta.out)
    modeller(nextclade.out[1])
    pangolin(manta.out)
    snpEff(vcf_for_fasta.out)
    simpleStats(manta.out.join(wgsMetrics.out))

    // Coinfection line
    coinfection_ivar(bwa.out, indexGenome.out, primers)
    freyja(coinfection_ivar.out[0], indexGenome.out)
    coinfection_varscan(coinfection_ivar.out[1])
    coinfection_analysis(coinfection_varscan.out)
}
