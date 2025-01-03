process copy_genome_and_primers {
  tag "Detecting type for sample:${sampleId}"
  container  = params.main_image
  // Dummy modules that pass genome from within container to the nextflow work dir
  // It is dummy hence there is no difference between illumina and nanopre
input:
  tuple val(sampleId), path(reads), val(QC_STATUS)
output:
  tuple val(sampleId), path("sars*"), path("primers.bed"), path("pairs.tsv"), env(REF_GENOME_ID), env(QC_exit), emit: all
  tuple val(sampleId), path("sars*"), env(QC_exit), emit: to_bwa
  tuple val(sampleId), path("primers.bed"), path("pairs.tsv"), emit: primers_and_pairs
  tuple val(sampleId), path("primers.bed"), emit: only_primers // required for nanopore second round
  tuple val(sampleId), path("genome.fasta"), emit: only_genome // indelqual module requires a variable not a tupple
  tuple val(sampleId), env(REF_GENOME_ID), emit: to_snpeff 
script:
"""

if [ ${QC_STATUS} == "nie" ]; then
  QC_exit="nie"
  REF_GENOME_ID="unk"
  touch pairs.tsv
  touch RSV_dummy.fasta
  touch RSV_dummy.fasta.amb
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
fi

"""

}

