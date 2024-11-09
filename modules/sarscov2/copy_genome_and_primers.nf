process copy_genome_and_primers {
  tag "Detecting type for sample:${sampleId}"
  container  = params.main_image
  // Dummy modules that pass genome from within container to the nextflow work dir
  // It is dummy hence there is no difference between illumina and nanopre
input:
  tuple val(sampleId), path(reads), val(QC_STATUS)
output:
  tuple val(sampleId), path("sars*"), path("primers.bed"), path("pairs.tsv"), env(REF_GENOME_ID), env(QC_exit), emit: all
  tuple val(sampleId), path("genome.fasta"), path("primers.bed"), env(REF_GENOME_ID), env(QC_exit), emit: all_nanopore
  tuple val(sampleId), path("sars*"), path("primers.bed"), env(QC_exit), emit: to_bwa
  tuple val(sampleId), path("primers.bed"), path("pairs.tsv"), emit: primers_and_pairs
  tuple val(sampleId), path("primers.bed"), emit: primers
  tuple val(sampleId), path("genome.fasta"), emit: only_genome // indelqual module requires a variable not a tupple
script:
"""

if [ ${QC_STATUS} == "nie" ]; then
  QC_exit="nie"
  REF_GENOME_ID="unk"
  touch pairs.tsv
  touch sars2_dummy.fasta
  touch sars2_dummy.fasta.amb
  touch primers.bed
  touch genome.fasta
else
  QC_exit="tak"
  cp /home/data/sarscov2/primers/${params.primers_id}/*bed primers.bed
  cp /home/data/sarscov2/primers/${params.primers_id}/*tsv pairs.tsv
  cp /home/data/sarscov2/genome/* .
  cp /home/data/sarscov2/genome/sarscov2.fasta genome.fasta
  REF_GENOME_ID=`head -1 genome.fasta | cut -d " " -f1 | tr -d ">"`
fi

"""

}

