// // // CHANGE THAT PRIOR TO JENKINS // // //
params.projectDir = ""
modules = "${params.projectDir}/modules" // Modules are part of the project_dir

// External databases with PREDEFINED structure
// When we use EXTERNAL database within a module we mount this path
// and the module itself will access the relevant database for a relevant species

params.external_databases_path=""

// All docker images used by this pipeline
// All modules now require explicit information which of these images they use
params.main_image = "nf_viral_main:1.1"
params.manta_image = "nf_viral_manta:1.1"
params.medaka_image = "ontresearch/medaka:sha447c70a639b8bcf17dc49b51e74dfcde6474837b-amd64"

// // // // END END END END END // // // // //



// Input section, only this section can be modified by shell wrapper for Pawel
// Other parameters should not be modifie to allow reporoducibility between samples
params.machine = '' // Can be set to either 'Illumina' or 'Nanopore'. Required
params.reads = "" // Must be provided by user
params.primers_id = "" // Organism specific, e.g V0 or V1 for RSV
params.adapters_id="TruSeq3-PE-2" // Can be user-provided, there is a default
params.species = "" // Required, name of a species expected to be analyzed with this pipeline can be either "SARS-CoV-2", "RSV" or "Influenza" 

// Output Directory
params.results_dir = "./results/"

// project_dir is not defined in sars pipeline explicite but is set by bash wrapper
// so we define it here, it can be used for version control
// Modules are just a subdirectory to project_dir and this cannot be changed

// All species-relevant variables, for now only expected genus for kraken2
// Furthermore if a user provides a wrong species the pipeline will not execute

if ( params.species  == 'SARS-CoV-2' ) {
genus="Betacoronavirus"
params.max_number_for_SV = 200000 
} else if (params.species  == 'Influenza') {
genus="Alphainfluenzavirus"
params.variant = "UNK"
params.max_number_for_SV = 10000
// Betainfluenzavirus for B/ kraken2 for now only undestands one genus

} else if (params.species  == 'RSV') {
genus="Orthopneumovirus"
params.max_number_for_SV = 100000

} else {
  println("Incorrect species, avalable options are : SARS-CoV-2, RSV or Influenza")
  System.exit(0)
}

// High-level parameters 
params.memory = 4048
params.threads = 5

// All sequencing platform specific parameters
// Some of them should be made common
if ( params.machine  == 'Illumina' ) {
params.min_number_of_reads = 1
params.expected_genus_value = 5
params.min_median_quality = 0
params.quality_initial = 5
params.length = 90
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
params.freyja_minq = 20
} else if (params.machine  == 'Nanopore') {
params.freyja_minq = 2
params.bed_offset=10 // for filtering
params.extra_bed_offset=10 // for filtering
params.min_mapq = 30 // for filtering
params.window_size = 50 // for filtering window size in which we equalize the coverage 
params.length = 0.49 // for filtering, nanopore min length is relative to the expected segment/amplikon length


params.medaka_model = "r941_min_sup_variant_g507" // Flow cell v9.4.1, for first round of medaka for the second round we use r941_min_sup_g507

if ( params.species  == 'SARS-CoV-2' ) {
params.medaka_chunk_len =  5000  // Flow cell v9.4.1
params.medaka_chunk_overlap = 4000 // Flow cell v9.4.1
} else if (params.species  == 'Influenza') {
// Betainfluenzavirus for B/ kraken2 for now only undestands one genus
params.medaka_chunk_len = 1000  // Flow cell v9.4.1
params.medaka_chunk_overlap = 500 // Flow cell v9.4.1

} else if (params.species  == 'RSV') {
params.medaka_chunk_len = 5000  // Flow cell v9.4.1
params.medaka_chunk_overlap = 4000 // Flow cell v9.4.1
}


params.min_number_of_reads = 1 // Stop the analysis if after mapping step bam has less than that number of reads

params.expected_genus_value = 5
params.min_median_quality = 0
params.quality_initial = 2 // We are extreamly lenient for nanopore 
params.max_depth = 600
params.min_cov = 50
params.mask = 50

params.quality_snp = 5 // We are extreamly lenient for nanopore

params.first_round_pval = 0.25 
params.second_round_pval = 0.05
params.pval = 0.05 // For varscan
params.lower_ambig = 0.45
params.upper_ambig = 0.55
params.quality_for_coverage = 1 // Parametr uzywany w modul lowCov, again we are extreamly lenient 

} else {
  println("Incorrect sequnecing platform, avalable options are : Illumina and Nanopore")
  System.exit(0)
}

// Modules section
// we load all the modules even if we do not plan to use all of them for a given species/platform combination

// More or less common modules
// In principle all modules could be made common if we write enough if-else statesment with them

include { kraken2_illumina } from "${modules}/common/kraken2.nf"
include { kraken2_nanopore } from "${modules}/common/kraken2.nf"

include { bwa } from "${modules}/common/bwa.nf"
include { minimap2 as minimap2_1 } from "${modules}/common/minimap2.nf"
include { minimap2 as minimap2_2 } from "${modules}/common/minimap2.nf"

include { dehumanization_illumina } from "${modules}/common/dehumanization.nf"
include { dehumanization_nanopore } from "${modules}/common/dehumanization.nf"

include { fastqc as fastqc_1 } from "${modules}/common/fastqc.nf"
include { run_fastqc_nanopore as run_fastqc_nanopore_1 } from "${modules}/common/fastqc_nanopore.nf"
include { fastqc as fastqc_2 } from "${modules}/common/fastqc.nf"

include { trimmomatic } from "${modules}/common/trimmomatic.nf"

include { filtering as filtering_one_segment } from "${modules}/sarscov2/filtering.nf" 

include { filtering_nanopore as filtering_one_segment_nanopore_1 } from "${modules}/sarscov2/filtering.nf"
include { filtering_nanopore as filtering_one_segment_nanopore_2 } from "${modules}/sarscov2/filtering.nf"

include { masking } from "${modules}/common/masking.nf"
include { masking_nanopore as masking_nanopore_strict_1 } from "${modules}/common/masking.nf"
include { masking_nanopore as masking_nanopore_strict_2 } from "${modules}/common/masking.nf"

include { masking_nanopore as masking_nanopore_overshot_1 } from "${modules}/common/masking.nf"
include { masking_nanopore as masking_nanopore_overshot_2 } from "${modules}/common/masking.nf"

include { merging } from "${modules}/sarscov2/merging.nf" 

include { merging_nanopore as merging_nanopore_1 } from "${modules}/common/merging_nanopore.nf"
include { merging_nanopore as merging_nanopore_2 } from "${modules}/common/merging_nanopore.nf"

include { picard_downsample_multisegment as picard_downsample } from "${modules}/common/picard.nf"

include { introduce_SV_with_manta } from "${modules}/common/manta.nf"

include { indelQual } from "${modules}/common/indelQual.nf"
include { picard_wgsMetrics } from "${modules}/common/picard_wgsMetrics.nf"
include { lowCov } from "${modules}/common/lowCov.nf"
include { varScan } from "${modules}/common/varscan.nf"
include { freeBayes } from "${modules}/common/freeBayes.nf"
include { lofreq } from "${modules}/common/lofreq.nf"
include { consensus_illumina } from "${modules}/common/consensus.nf"
include { consensus_nanopore } from "${modules}/common/consensus.nf"

// The output for different species/platforms has different channels, hence we need to define
// 4 different versions of aggregator module
include { json_aggregator_sars_illumina } from "${modules}/common/json_aggregator.nf"
include { json_aggregator_sars_nanopore } from "${modules}/common/json_aggregator.nf"
include { json_aggregator_rsv_illumina } from "${modules}/common/json_aggregator.nf"
include { json_aggregator_rsv_nanopore } from "${modules}/common/json_aggregator.nf"
include { json_aggregator_influenza_illumina } from "${modules}/common/json_aggregator.nf"
include { json_aggregator_influenza_nanopore } from "${modules}/common/json_aggregator.nf"

// vcf_for_fasta, snpEff and nextclade should be "common" module ? check if they work for influenza 
include { vcf_for_fasta } from "${modules}/sarscov2/vcf_for_fasta.nf"
include { snpEff_illumina } from "${modules}/sarscov2/snpEff.nf"
include { snpEff_nanopore } from "${modules}/sarscov2/snpEff.nf"

include { nextclade as nextclade_noninfluenza } from "${modules}/sarscov2/nextclade.nf"

include { modeller } from "${modules}/sarscov2/modeller.nf"

// SARS specific modules
include { copy_genome_and_primers } from "${modules}/sarscov2/copy_genome_and_primers.nf"
include { pangolin } from "${modules}/sarscov2/pangolin.nf"

// // Coinfection line for SARS
include { freyja as freyja_sars } from "${modules}/sarscov2/freyja.nf"
include { coinfection_genome_masking_illumina } from "${modules}/sarscov2/coinfection_ivar.nf"
include { coinfection_genome_masking_nanopore } from "${modules}/sarscov2/coinfection_ivar.nf"

include { coinfection_varscan as coinfection_varscan_sars } from "${modules}/sarscov2/coinfection_varscan.nf"
include { coinfection_analysis_illumina as coinfection_analysis_illumina_sars } from "${modules}/sarscov2/coinfection_analysis.nf"
include { coinfection_analysis_nanopore as coinfection_analysis_nanopore_sars } from "${modules}/sarscov2/coinfection_analysis.nf"

// INFL-specific modules 
include { detect_subtype_illumina as detect_subtype_influenza_illumina } from "${modules}/infl/detect_subtype.nf"
include { detect_subtype_nanopore as detect_subtype_influenza_nanopore } from "${modules}/infl/detect_subtype_nanopore.nf"


include { reassortment as reassortment_influenza } from "${modules}/infl/reassortment.nf"
include { filtering as filtering_influenza_illumina} from "${modules}/infl/filtering.nf"

include { filtering_nanopore as filtering_influenza_nanopore_1 } from "${modules}/infl/filtering.nf"
include { filtering_nanopore as filtering_influenza_nanopore_2 } from "${modules}/infl/filtering.nf"

include { sort_and_index as sort_and_index_influenza_illumina } from "${modules}/infl/sort_and_index.nf"


include { nextclade as nextclade_influenza } from "${modules}/infl/nextclade.nf"
include { nextalign as nextalign_influenza } from "${modules}/infl/nextalign.nf"
include { resistance as resistance_influenza } from "${modules}/infl/resistance.nf"

// RSV-specific modules
include { detect_type_illumina as detect_type_rsv_illumina } from "${modules}/rsv/detect_type.nf"
include { detect_type_nanopore as detect_type_rsv_nanopore } from "${modules}/rsv/detect_type.nf"

include { medaka_first_round as medaka_1 } from "${modules}/common/medaka.nf"
include { medaka_second_round as medaka_2 } from "${modules}/common/medaka.nf"

include { medaka_varscan_integration_first_round } from "${modules}/common/integrate_medaka_and_varscan.nf"
include { medaka_varscan_integration_second_round } from "${modules}/common/integrate_medaka_and_varscan.nf"

include { filter_out_non_SNPs as filter_out_non_SNPs_1 } from "${modules}/common/filter_out_non_SNPs.nf"

include { make_genome_from_vcf as make_genome_from_vcf_1 } from "${modules}/common/make_genome_from_vcf.nf"
include { make_genome_from_vcf as make_genome_from_vcf_2 } from "${modules}/common/make_genome_from_vcf.nf"

include { varScan as varScan_1 } from "${modules}/common/varscan.nf"
include { varScan as varScan_2 } from "${modules}/common/varscan.nf"

include { substitute_ref_genome } from "${modules}/sarscov2/substitute_ref.nf"

// Script execution path
ExecutionDir = new File('.').absolutePath // all output in json should be relative to this path

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
        detect_type_out = copy_genome_and_primers(final_reads_and_final_qc)
    } else if (params.species  == 'Influenza') {
        detect_subtype_out = detect_subtype_influenza_illumina(final_reads_and_final_qc)
        detect_type_out =  reassortment_influenza(detect_subtype_out.segments_scores)
    } else if (params.species  == 'RSV') {
        detect_type_out = detect_type_rsv_illumina(final_reads_and_final_qc)
    }
    reads_and_genome = trimmomatic_out.proper_reads.join(detect_type_out.to_bwa, by:0)
    // Mapping
    bwa_out = bwa(reads_and_genome)

    // Dehumanization
    dehumanization_out = dehumanization_illumina(bwa_out.only_bam.join(trimmomatic_out.proper_reads, by:0))

    // coinfection analysis for SARS ONLY !
    if ( params.species  == 'SARS-CoV-2' ) {
        coinfection_ivar_sars_out = coinfection_genome_masking_illumina(bwa_out.to_coinfection)
        freyja_out =  freyja_sars(coinfection_ivar_sars_out.to_freyja)
        coinfection_varscan_out = coinfection_varscan_sars(coinfection_ivar_sars_out.to_custom_analysis)
        coinfection_analysis_sars_out = coinfection_analysis_illumina_sars(coinfection_varscan_out)
        coinfection_json = freyja_out.json.join(coinfection_analysis_sars_out.json)
    }

    // Filtering script is different for one-segment and multiple-segments organisms
    if ( params.species  == 'SARS-CoV-2' || params.species  == 'RSV' ) {
        initial_bam_and_primers = bwa_out.only_bam.join(detect_type_out.primers_and_pairs, by:0)
        filtering_out = filtering_one_segment(initial_bam_and_primers)
        masking_out = masking(filtering_out.one_amplicon_primers_and_QC)
        merging_out = merging(filtering_out.two_amplicon_only.join(masking_out, by:0))
        pre_final_bam_and_genome = merging_out.join(detect_type_out.only_genome, by:0)
    } else if (params.species  == 'Influenza') {
        initial_bam_and_primers = bwa_out.bam_and_genome.join(detect_type_out.primers_and_pairs, by:0)
        filtering_out = filtering_influenza_illumina(initial_bam_and_primers)
        masking_out = masking(filtering_out.one_amplicon_primers_and_QC)
        sort_and_index_out = sort_and_index_influenza_illumina(masking_out)
        pre_final_bam_and_genome = sort_and_index_out.join(detect_type_out.only_genome, by:0)
    }

    // Predicting sample's genome is identical
    indelQual_out = indelQual(pre_final_bam_and_genome)
    lowCov_out = lowCov(indelQual_out.bam_genome_and_qc)
    varScan_out = varScan(indelQual_out.bam_genome_and_qc)
    freebayes_out = freeBayes(indelQual_out.bam_genome_and_qc)
    lofreq_out = lofreq(indelQual_out.bam_genome_and_qc)
    wgsMetrics_out = picard_wgsMetrics(indelQual_out.bam_genome_and_qc.join(filtering_out.json))
    all_sub_fastas = lowCov_out.fasta.join(varScan_out.fasta)
    all_sub_fastas = all_sub_fastas.join(freebayes_out)
    all_sub_fastas = all_sub_fastas.join(lofreq_out)
    consensus_out = consensus_illumina(all_sub_fastas)

    // Predicting SV with manta
    picard_downsample_out = picard_downsample(bwa_out.bam_and_genome)
    final_genome_out = introduce_SV_with_manta(picard_downsample_out.to_manta.join(consensus_out.multiple_fastas, by:0))


  } else if (params.machine == 'Nanopore') {
        Channel
            .fromPath(params.reads)
            .map {it -> tuple(it.getName().split("\\.")[0..<(-2)].join('_'), it)}
            .set {reads}
        fastqc_initial_out = run_fastqc_nanopore_1(reads, "pre-filtering")

        reads_and_qc = reads.join(fastqc_initial_out.qcstatus)
        kraken2_out = kraken2_nanopore(reads_and_qc, genus)
        final_reads_and_final_qc = reads.join(kraken2_out.qcstatus_only, by:0)
        

        // Get initial reference genome and primers
        if ( params.species  == 'SARS-CoV-2' ) {
          // For SARS2 we only need to copy data
          detect_type_out = copy_genome_and_primers(final_reads_and_final_qc)
        } else if (params.species  == 'Influenza') {
          detect_subtype_out = detect_subtype_influenza_nanopore(final_reads_and_final_qc)
          detect_type_out =  reassortment_influenza(detect_subtype_out.segments_scores)
        } else if (params.species  == 'RSV') {
          detect_type_out = detect_type_rsv_nanopore(final_reads_and_final_qc)
        }
  
        reads_and_genome_and_primers = reads.join(detect_type_out.all_nanopore, by:0)
 
        minimap2_1_out = minimap2_1(reads_and_genome_and_primers)
      
        if ( params.species  == 'SARS-CoV-2' ) {
          // kanal minimap2_1 jest itenconalnie
          coinfection_genome_masking_out = coinfection_genome_masking_nanopore(minimap2_1_out.bam_and_genome_and_primers)
          freyja_out =  freyja_sars(coinfection_genome_masking_out.to_freyja)
          coinfection_varscan_out = coinfection_varscan_sars(coinfection_genome_masking_out.to_custom_analysis)
          coinfection_analysis_sars_out = coinfection_analysis_nanopore_sars(coinfection_varscan_out)
          coinfection_json = freyja_out.json.join(coinfection_analysis_sars_out.json)
        }
 
        if ( params.species  == 'SARS-CoV-2' ||  params.species  == 'RSV') {
          filteing_1_out = filtering_one_segment_nanopore_1(minimap2_1_out.bam_and_genome_and_primers)
          normal_masking_1_out = masking_nanopore_strict_1(filteing_1_out.to_normal_masking, params.bed_offset, 0)
          overshot_masking_1_out = masking_nanopore_overshot_1(filteing_1_out.to_overshot_masking, params.bed_offset + params.extra_bed_offset, 0)
          merging_1_out = merging_nanopore_1(normal_masking_1_out.bam_and_genome.join(overshot_masking_1_out.bam_only))
          to_medaka_1 = merging_1_out.to_medaka
        } else if (params.species  == 'Influenza') {
          filtering_1_out = filtering_influenza_nanopore_1(minimap2_1_out.bam_and_genome_and_primers)
          normal_masking_1_out = masking_nanopore_strict_1(filtering_1_out.to_normal_masking, params.bed_offset, 0)
          to_medaka_1 = normal_masking_1_out.bam_and_genome
        }

        // First round of medaka
        // Rough approximation of a final genome
        medaka_1_out = medaka_1(to_medaka_1)
        varScan_1_out = varScan_1(to_medaka_1)

        medaka_varscan_integration_1_out = medaka_varscan_integration_first_round(medaka_1_out.vcf.join(varScan_1_out.pre_vcf))
        filter_out_non_SNPs_1_out = filter_out_non_SNPs_1(medaka_varscan_integration_1_out.vcf)
        novel_genome =  make_genome_from_vcf_1(medaka_varscan_integration_1_out.reference_genome.join(filter_out_non_SNPs_1_out.vcf))
 
        // Second round of medaka using genome obtained in round 1
        // at the end of this round we need need yo resote the origunal genome as reference (for snpeff)

        reads_and_updated_genome_and_primers = reads.join(novel_genome.only_fasta, by:0)
        reads_and_updated_genome_and_primers = reads_and_updated_genome_and_primers.join(detect_type_out.primers,  by:0)
        reads_and_updated_genome_and_primers = reads_and_updated_genome_and_primers.join(novel_genome.only_QC, by: 0)

        minimap2_2_out = minimap2_2(reads_and_updated_genome_and_primers)

        dehumanization_out = dehumanization_nanopore(minimap2_2_out.bam_and_qc.join(reads, by:0))
       
        if ( params.species  == 'SARS-CoV-2' ||  params.species  == 'RSV') {
          filtering_2_out = filtering_one_segment_nanopore_2(minimap2_2_out.bam_and_genome_and_primers)
          normal_masking_2_out = masking_nanopore_strict_2(filtering_2_out.to_normal_masking, params.bed_offset, 0)
          overshot_masking_2_out = masking_nanopore_overshot_2(filtering_2_out.to_overshot_masking, params.bed_offset + params.extra_bed_offset, 0)
          merging_2_out = merging_nanopore_2(normal_masking_2_out.bam_and_genome.join(overshot_masking_2_out.bam_only))
          to_medaka_2 = merging_2_out.to_medaka
        } else if (params.species  == 'Influenza') {
          filtering_2_out = filtering_influenza_nanopore_2(minimap2_2_out.bam_and_genome_and_primers)
          normal_masking_2_out = masking_nanopore_strict_2(filtering_2_out.to_normal_masking, params.bed_offset, 0)
          to_medaka_2 = normal_masking_2_out.bam_and_genome
        }

        wgsMetrics_out = picard_wgsMetrics(to_medaka_2.join(filtering_2_out.json))
        medaka_2_out = medaka_2(to_medaka_2)
        varScan_2_out = varScan_2(to_medaka_2)

        medaka_varscan_integration_2_out = medaka_varscan_integration_second_round(medaka_2_out.vcf.join(varScan_2_out.pre_vcf))
        novel_genome_2_out = make_genome_from_vcf_2(medaka_varscan_integration_2_out.vcf)
        lowCov_out = lowCov(to_medaka_2)
        to_final_genome = lowCov_out.fasta.join(novel_genome_2_out.fasta_and_QC)
        to_final_genome = to_final_genome.join(medaka_varscan_integration_2_out.reference_genome)
        prefinal_genome_out = consensus_nanopore(to_final_genome)  
        final_genome_out = substitute_ref_genome(prefinal_genome_out.fasta_refgenome_and_qc.join(detect_type_out.only_genome))
 
  }
  // Post FASTA generation modules mostly common for nanopore and illumina

  if ( params.species  == 'SARS-CoV-2' || params.species  == 'RSV' ) {
      nextclade_out = nextclade_noninfluenza(final_genome_out.fasta_refgenome_and_qc)
  } else if (params.species  == 'Influenza') {
      // manta_out.fasta_refgenome_and_qc.join(detect_subtype_illumina_out.subtype_id, by:0)
      final_genome_and_influenza_subtype = final_genome_out.fasta_refgenome_and_qc.join(detect_subtype_out.subtype_id, by:0)
      nextclade_out = nextclade_influenza(final_genome_and_influenza_subtype)
      nextalign_out = nextalign_influenza(final_genome_and_influenza_subtype)
      resistance_out = resistance_influenza(nextalign_out)
  }
  modeller_out = modeller(nextclade_out.to_modeller)
  // Pangolin only for SARS, module is species-aware
  pangolin_out = pangolin(final_genome_out.fasta_refgenome_and_qc)

  // final vcf + snpEFF, snpEFF is species-aware
  vcf_for_fasta_out = vcf_for_fasta(final_genome_out.fasta_refgenome_and_qc)
   
  // illumina/nanopore specific channel 
  if(params.machine == 'Illumina') {
    snpEff_out = snpEff_illumina(vcf_for_fasta_out.vcf.join(indelQual_out.bam_and_qc, by:0))
  }
  else if (params.machine == 'Nanopore') {
    snpEff_out = snpEff_nanopore(vcf_for_fasta_out.vcf.join(minimap2_2_out.bam_and_qc, by:0))  
  }


  // JSON OUTPUT SECTION
  for_json_aggregator = fastqc_initial_out.json.join(kraken2_out.json)

  if(params.machine == 'Illumina') {
    for_json_aggregator = for_json_aggregator.join(fastqc_filtered_out.json) // only fo illumina for nanopore dummy process is used to prodeuce json channel
  }

  if ( params.species  == 'Influenza' ) {
    for_json_aggregator = for_json_aggregator.join(detect_type_out.json)
  }

  if(params.machine == 'Illumina') {
    for_json_aggregator = for_json_aggregator.join(bwa_out.json)
  } else if (params.machine == 'Nanopore') {
    for_json_aggregator = for_json_aggregator.join(minimap2_1_out.json)
  }

  if ( params.species  == 'SARS-CoV-2' ) {
    for_json_aggregator = for_json_aggregator.join(coinfection_json)
  }

  for_json_aggregator = for_json_aggregator.join(dehumanization_out.json) 
  for_json_aggregator = for_json_aggregator.join(wgsMetrics_out.json)
 

  if(params.machine == 'Illumina') {
    for_json_aggregator = for_json_aggregator.join(final_genome_out.json)
  } else if (params.machine == 'Nanopore') {
    for_json_aggregator = for_json_aggregator.join(prefinal_genome_out.json) // tylko nanopore
  }
  
  for_json_aggregator = for_json_aggregator.join(pangolin_out.json)
  for_json_aggregator = for_json_aggregator.join(nextclade_out.json)
  for_json_aggregator = for_json_aggregator.join(snpEff_out)
  for_json_aggregator = for_json_aggregator.join(modeller_out.json)
  
  if ( params.species  == 'Influenza' ) {
    for_json_aggregator = for_json_aggregator.join(resistance_out.json)
  }


  if(params.species == 'SARS-CoV-2') {
    if(params.machine == 'Illumina') {
      json_aggregator_sars_illumina(for_json_aggregator, ExecutionDir)
    } else if (params.machine == 'Nanopore') {
      json_aggregator_sars_nanopore(for_json_aggregator, ExecutionDir)
    }
  
  } else if (params.species == 'RSV') {
    if(params.machine == 'Illumina') {
      json_aggregator_rsv_illumina(for_json_aggregator, ExecutionDir)
    } else if (params.machine == 'Nanopore') {
      json_aggregator_rsv_nanopore(for_json_aggregator, ExecutionDir)
    }
  } else if (params.species == 'Influenza') {
    if(params.machine == 'Illumina') {
      json_aggregator_influenza_illumina(for_json_aggregator, ExecutionDir)
    } else if (params.machine == 'Nanopore') {
      json_aggregator_influenza_nanopore(for_json_aggregator, ExecutionDir)
    }

  }
}
