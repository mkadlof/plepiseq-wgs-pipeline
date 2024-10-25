// Input section, only this section can be modified by shell wrapper for Pawel
// Other parameters should not be modifie to allow reporoducibility between samples
params.machine = '' // Can be set to either 'Illumina' or 'Nanopore'. Required
params.reads = "" // Must be provided by user
params.primers_id = "" // Either V0 or V1
params.adapters_id="TruSeq3-PE-2" // To podaje user ale jest default
params.species = "" // Required field name of a species expected to be analyzed with this pipeline can be either "SARS-CoV-2", "RSV" or "Influenza" 

// Output Dir by default "results" directory 
params.results_dir = "./results"

// projectDir is not defined in sars pipeline in nf explicite but is set by bash wrapper
// Modules are just a subdirectory to projectDir and this cannot be changed

params.projectDir = "/home/michall/git/nf_illumina_sars_ml"
modules = "${params.projectDir}/modules" // Modules are part of the projectDir

// External databases with PREDEFINED structure
// When we use EXTERNAL database with a module we mount this path and the module itself will access the relevant database for relevant species
params.external_databases_path="/home/jenkins/workspace/nf_illumina_sars/external_databases"

// All species-relevant variables, for now only expected genus for kraken2
if ( params.species  == 'SARS-CoV-2' ) {
genus="Betacoronavirus"

} else if (params.species  == 'Influenza') {
genus="Alphainfluenzavirus"
// Betainfluenzavirus for B/ kraken2 for now only undestands one genus

} else if (params.species  == 'RSV') {
genus="Orthopneumovirus"

} else {
  println("Incorrect species, avalable options are : SARS-CoV-2, RSV or Influenza")
  System.exit(0)
}

// All docker images used by this pipeline
// All modules require explicit information which of this images they use
params.main_image = "nf_viral_main:1.0"
params.manta_image = "nf_viral_manta:1.0"
params.medaka_image = "ontresearch/medaka:sha447c70a639b8bcf17dc49b51e74dfcde6474837b-amd64"

// High-level parameters 
params.memory = 4048
params.threads = 5

// All sequencing platform specific parameters
// Some of them should be made common
if ( params.machine  == 'Illumina' ) {
params.min_number_of_reads = 0
params.expected_genus_value = 5
params.min_median_quality = 0
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
params.quality_for_coverage = 10 // Parametr uzywany w modul lowCov
} else if (params.machine  == 'Nanopore') {
params.model_medaka = "r941_min_hac_g507" // Flow cell v9.4.1
params.min_number_of_reads = 0
params.expected_genus_value = 5
params.min_median_quality = 0
params.quality_initial = 2 // We are extreamly liberal for nanopore 
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
params.quality_for_coverage = 10 // Parametr uzywany w modul lowCov
} else {
  println("Incorrect sequnecing platform, avalable options are : Illumina and Nanopore")
  System.exit(0)
}

// Modules section, we can load all module even if we do not plan to use all of them for a given species/platform combination

// More or less common modules
// In principle all modules could be made common if we write enough if-else statesment with them
include { kraken2_illumina } from "${modules}/common/kraken2.nf"
include { kraken2_nanopore } from "${modules}/common/kraken2.nf"

include { bwa } from "${modules}/common/bwa.nf"
include { minimap2 } from "${modules}/common/minimap2.nf"
include { dehumanization_illumina } from "${modules}/common/dehumanization.nf"
include { dehumanization_nanopore } from "${modules}/common/dehumanization.nf"

include { fastqc as fastqc_1 } from "${modules}/common/fastqc.nf"
include { run_fastqc_nanopore as run_fastqc_nanopore_1 } from "${modules}/common/fastqc_nanopore.nf"
include { fastqc as fastqc_2 } from "${modules}/common/fastqc.nf"

include { trimmomatic } from "${modules}/common/trimmomatic.nf"

include { filtering as filtering_one_segment } from "${modules}/sarscov2/filtering.nf" 

include { masking } from "${modules}/common/masking.nf"
include { merging } from "${modules}/sarscov2/merging.nf" 
include { picard_downsample } from "${modules}/common/picard.nf"

include { introduce_SV_with_manta } from "${modules}/common/manta.nf"

include { indelQual } from "${modules}/common/indelQual.nf"
include { picard_wgsMetrics } from "${modules}/common/wgsMetrics.nf"
include { lowCov } from "${modules}/common/lowCov.nf"
include { varScan } from "${modules}/common/varscan.nf"
include { freeBayes } from "${modules}/common/freeBayes.nf"
include { lofreq } from "${modules}/common/lofreq.nf"
include { consensus } from "${modules}/common/consensus.nf"

// vcf_for_fasta, snpEff and nextclade should be "common" module ? check if they work for influenza 
include { vcf_for_fasta } from "${modules}/sarscov2/vcf_for_fasta.nf"
include { snpEff } from "${modules}/rsv/snpEff.nf"
include { nextclade } from "${modules}/sarscov2/nextclade.nf"

// include { simpleStats } from "${params.modules}/sarscov2/simpleStats.nf"
// include { modeller } from "${params.modules}/sarscov2/modeller.nf"
// include { json_aggregator } from "${params.modules}/common/json_aggregator.nf"

// SARS specific modules
include { copy_genome_and_primers } from "${modules}/sarscov2/copy_genome_and_primers.nf"
include { pangolin } from "${modules}/sarscov2/pangolin.nf"

// // Coinfection line for SARS
include { freyja } from "${params.modules}/sarscov2/freyja.nf"
include { coinfection_ivar } from "${params.modules}/sarscov2/coinfection_ivar.nf"
include { coinfection_varscan } from "${params.modules}/sarscov2/coinfection_varscan.nf"
include { coinfection_analysis } from "${params.modules}/sarscov2/coinfection_analysis.nf"

// INFL-specific modules 
include { detect_subtype as detect_subtype_influenza_illuumina } from "${modules}/infl/detect_subtype.nf"
include { reassortment as reassortment_influenza_illumina } from "${modules}/infl/reassortment.nf"
include { filtering as filtering_multiple_segments} from "${modules}/infl/filtering.nf"
include { sort_and_index as sort_and_index_influenza_illumina } from "${modules}/infl/sort_and_index.nf"
include { nextclade as nextclade_influenza } from "${modules}/infl/nextclade.nf"
include { nextalign } from "${modules}/infl/nextalign.nf"
include { resistance } from "${modules}/infl/resistance.nf"

// RSV-specific modules
include { detect_type_rsv_illumina } from "${modules}/rsv/detect_type.nf"
include { detect_type_rsv_nanopore } from "${modules}/rsv/detect_type.nf"

// Main workflow
workflow{

if(params.machine == 'Illumina') {
  Channel
    .fromFilePairs(params.reads)
    .set {reads}
  // Initail fastqc
  fastqc_initial_out = fastqc_1(reads, "pre-filtering")

  // Running kraken2 prediction
  reads_and_qc = reads.join(fastqc_initial_out.qcstatus)
  kraken2_out = kraken2_illumina(reads_and_qc, genus) 
  trimmomatic_out = trimmomatic(reads.join(kraken2_out.qcstatus_only, by:0))
  fastqc_filtered_out = fastqc_2(trimmomatic_out.proper_reads_and_qc, "post-filtering")
  
  final_reads_and_final_qc = trimmomatic_out.proper_reads.join(fastqc_filtered_out.qcstatus, by:0)
  
  // For all three species selection of reference genome/primers is different 
  if ( params.species  == 'SARS-CoV-2' ) {
    detect_type_illumina_out = copy_genome_and_primers(final_reads_and_final_qc)
    reads_and_genome = trimmomatic_out.proper_reads.join(detect_type_illumina_out.to_bwa, by:0)
  
  } else if (params.species  == 'Influenza') {
    // TBD
  } else if (params.species  == 'RSV') {
    detect_type_illumina_out = detect_type_rsv_illumina(final_reads_and_final_qc)
    reads_and_genome = trimmomatic_out.proper_reads.join(detect_type_illumina_out.to_bwa, by:0)
  }
  
 
  bwa_out = bwa(reads_and_genome)

  // Dehumanization
  dehumanization_illumina_out = dehumanization_illumina(bwa_out.only_bam.join(trimmomatic_out.proper_reads, by:0))

  initial_bam_and_primers = bwa_out.only_bam.join(detect_type_illumina_out.primers_and_pairs, by:0)
  
  filtering_out = filtering(initial_bam_and_primers)
  masking_out = masking(filtering_out.one_amplicon_primers_and_QC)
  merging_out = merging(filtering_out.two_amplicon_only.join(masking_out, by:0)) 

  pre_final_bam_and_genome = merging_out.join(detect_type_illumina_out.only_genome, by:0)
  indelQual_out = indelQual(pre_final_bam_and_genome) 
 
  final_bam_and_genome = indelQual_out.bam_and_qc.join(detect_type_illumina_out.only_genome, by:0)   
  
  freebayes_out = freeBayes(final_bam_and_genome)
  lofreq_out = lofreq(final_bam_and_genome)
  lowCov_out = lowCov(final_bam_and_genome) 
  varScan_out = varScan(final_bam_and_genome)
  wgsMetrics_out = picard_wgsMetrics(final_bam_and_genome) 

  all_sub_fastas = lowCov_out.fasta.join(varScan_out)
  all_sub_fastas = all_sub_fastas.join(freebayes_out)
  all_sub_fastas = all_sub_fastas.join(lofreq_out)
  
  consensus_out = consensus(all_sub_fastas) 
  genome_sequence_and_ref = consensus_out.single_fasta.join(detect_type_illumina_out.only_genome, by:0)
  vcf_for_fasta_out = vcf_for_fasta(genome_sequence_and_ref)

  snpEff_out = snpEff(vcf_for_fasta_out.vcf.join(indelQual_out.bam_genome_and_qc, by:0))
  
  // Predicting SV with manta
  picard_downsample_out = picard_downsample(bwa_out.bam_and_genome)
  manta_out = introduce_SV_with_manta(picard_downsample_out.to_manta.join(consensus_out.single_fasta, by:0))
  nextclade_out = nextclade(manta_out.fasta_refgenome_and_qc)
 
} else if (params.machine == 'Nanopore') {
  Channel
  .fromPath(params.reads)
  .map {it -> tuple(it.getName().split("\\.")[0], it)}
  .set {reads}
 
  fastqc_initial_out = run_fastqc_nanopore_1(reads, "pre-filtering")
 
  reads_and_qc = reads.join(fastqc_initial_out.qcstatus)
  kraken2_out = kraken2_nanopore(reads_and_qc, "Orthopneumovirus")
  final_reads_and_final_qc = reads.join(kraken2_out.qcstatus_only, by:0)
  detect_type_nanopore_out = detect_type_rsv_nanopore(final_reads_and_final_qc)
  reads_and_genome = reads.join(detect_type_nanopore_out.to_minimap2, by:0)
  minimap2_out = minimap2(reads_and_genome)

  // Dehumanization
  dehumanization_nanopore_out = dehumanization_nanopore(minimap2_out.join(reads, by:0))

}

    // Channels
    // ref_genome = Channel.fromFilePairs("${projectDir}/data/sarscov2/genome/*sarscov2.fasta", size: -1).collect().sort()
    // ref_genome_with_index = Channel.fromFilePairs("${projectDir}/data/sarscov2/genome/*sarscov2.fasta{,.amb,.ann,.bwt,.fai,.pac,.sa}", size: -1).collect().sort()
    // primers_and_pairs = Channel.fromFilePairs("${projectDir}/data/sarscov2/primers/${params.primers_id}/*{nCoV-2019.scheme.bed,pairs.tsv}").collect()
    // coinfections = Channel.fromPath("${projectDir}/data/sarscov2/coinfections/").first()
    // vcf_template = Channel.fromPath("${projectDir}/data/sarscov2/vcf_template/").first()
    // modeller_data = Channel.fromPath("${projectDir}/data/sarscov2/modeller/").first()

    // The following construct is needed to add the sampleId to the channels
    // with constant data (like reference genome or primers). This is done to
    // generalize cases where the reference genome may not be identical for
    // each sample (as is the case in the influenza pipeline, for example).
    // ref_genome = ref_genome
    // .combine(reads).map { _, genome, sampleId, reads ->
    //     return [sampleId, genome]
    //}
    // ref_genome_with_index = ref_genome_with_index
    // .combine(reads).map { _, genome_and_index_files, sampleId, reads ->
    //     return [sampleId, genome_and_index_files]
    // }
    // primers_and_pairs = primers_and_pairs
    // .combine(reads).map { _, primers_and_pairs, sampleId, reads ->
    //    return [sampleId, primers_and_pairs]
    // }

    // primers = primers_and_pairs.map{ sampleId, files -> [sampleId, files[0]] }

    // The following two variables are used exclusively to include pipeline version information in the resulting output.json file.
    // def repo_path = workflow.projectDir
    // version = Channel.value("git -C ${repo_path} rev-parse HEAD".execute().text.trim().substring(0, 7))
    // pathogen = Channel.value('sars2')

    // Processes
    // fastqc_1(reads, "initialfastq")
    // c1 = reads.join(fastqc_1.out[2])
    // kraken2(c1)
    // trimmomatic(reads, adapters)
    // fastqc_2(trimmomatic.out[0], "aftertrimmomatic")
    // bwa(trimmomatic.out[0].join(ref_genome_with_index))
    // c2 = bwa.out.join(trimmomatic.out[1])
    // dehumanization(c2)
    // c3 = bwa.out.join(primers)
    // filtering(c3)
    // c4 = filtering.out[0].join(primers_and_pairs)
    // masking(c4)
    // c5 = filtering.out[1].join(masking.out)
    // merging(c5)
    // picard(bwa.out)
    // indelQual(merging.out.join(ref_genome))
    // lowCov(indelQual.out.join(ref_genome))
    // varScan(indelQual.out.join(ref_genome))
    // freeBayes(indelQual.out.join(ref_genome))
    // lofreq(indelQual.out.join(ref_genome_with_index))
    // wgsMetrics(indelQual.out.join(ref_genome))
    // c6 = lowCov.out[1].join(varScan.out).join(freeBayes.out).join(lofreq.out)
    // consensus(c6)
    // c7 = consensus.out[0].join(ref_genome)
    // vcf_for_fasta(c7, vcf_template)
    // manta(picard.out.join(consensus.out[0]))
    // nextclade(manta.out)
    // modeller(nextclade.out[1], modeller_data)
    // pangolin(manta.out)
    // c8 = vcf_for_fasta.out.join(indelQual.out).join(ref_genome)
    // snpEff(c8)
    // c9 = manta.out.join(wgsMetrics.out[0]).join(primers)
    // simpleStats(c9)

    // Coinfection line
    // c10 = bwa.out.join(ref_genome).join(primers)
    // coinfection_ivar(c10)
    // c11 = coinfection_ivar.out[0].join(ref_genome)
    // freyja(c11)
    // coinfection_varscan(coinfection_ivar.out[1])
    // coinfection_analysis(coinfection_varscan.out, coinfections)


    // json_aggregator(pathogen, version, wgsMetrics.out[1])
}
