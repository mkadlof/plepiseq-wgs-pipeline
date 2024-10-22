process detect_type_illumina {
  tag "Detecting type for sample:${sampleId}"
  maxForks 5
// Determine if we are dealing with RSV A or B based on score ratio after mapping to two genomes
input:
  tuple val(sampleId), path(reads), val(QC_STATUS)
output:
  tuple val(sampleId), path("RSV*"), path("primers.bed"), path("pairs.tsv"), env(TYPE), env(REF_GENOME_ID), env(QC_exit), emit: all
  tuple val(sampleId), path("RSV*"), env(QC_exit), emit: to_bwa
  tuple val(sampleId), path("primers.bed"), path("pairs.tsv"), emit: primers_and_pairs
  tuple val(sampleId), path("genome.fasta"), emit: only_genome // indelqual module requires a variable not a tupple
  tuple val(sampleId), env(TYPE), emit: json
script:
"""
# genome jest w kontenerze

if [ ${QC_STATUS} == "nie" ]; then
  QC_exit="nie"
  TYPE="unk"
  REF_GENOME_ID="unk"
  touch pairs.tsv
  touch RSV_dummy.fasta
  touch RSV_dummy.fasta.amb
  touch primers.bed
  touch genome.fasta
else
  REFERENCE_GENOME_FASTA="/home/data/rsv/genome/RSV/RSV.fasta"
  bwa mem -t ${params.threads} \
          -T 30 \
          "\${REFERENCE_GENOME_FASTA}" \
          ${reads[0]} \
          ${reads[1]}  | \
          samtools view -@ ${params.threads} -Sb -o type_determination.bam -f 3 -F 2048 -

  ILE_A=`samtools view type_determination.bam | grep "hRSV/A/England/397/2017" | wc -l `
  ILE_B=`samtools view type_determination.bam | grep "hRSV/B/England/RE20000104/2020" | wc -l `

  if [[ \${ILE_A} -lt 1000 && \${ILE_B} -lt 1000 ]]; then
      QC_exit="nie"
      TYPE="unk"
      REF_GENOME_ID="unk"
      touch pairs.tsv
      touch RSV_dummy.fasta
      touch RSV_dummy.fasta.amb
      touch primers.bed
      touch genome.fasta
  elif [[ \${ILE_A} -gt 1000 || \${ILE_B} -gt 1000 ]]; then
    QC_exit="tak"
    if  [ \${ILE_A} -gt \${ILE_B} ]; then
      TYPE="A"
      cp /home/data/rsv/primers/A/${params.primers_id}/*bed primers.bed
      cp /home/data/rsv/primers/A/${params.primers_id}/*tsv pairs.tsv
      cp /home/data/rsv/genome/RSV_A/RSV* .
      cp /home/data/rsv/genome/RSV_A/RSV_A.fasta genome.fasta
    else
      TYPE="B"
      cp /home/data/rsv/primers/B/${params.primers_id}/*bed primers.bed
      cp /home/data/rsv/primers/B/${params.primers_id}/*tsv pairs.tsv
      cp /home/data/rsv/genome/RSV_B/RSV* .
      cp /home/data/rsv/genome/RSV_B/RSV_B.fasta genome.fasta
    
    fi
  REF_GENOME_ID=`head -1 genome.fasta | cut -d " " -f1 | tr -d ">"`
  fi
fi

"""

}

process detect_type_nanopore {
  tag "Detecting type for sample:${sampleId}"
  maxForks 5
// Determine if we are dealing with RSV A or B based on score ratio after mapping to two genomes
input:
  tuple val(sampleId), path(reads), val(QC_STATUS)
output:
  tuple val(sampleId), path("RSV*"), path("primers.bed"), path("pairs.tsv"), env(TYPE), env(REF_GENOME_ID), env(QC_exit), emit: all
  tuple val(sampleId), path("RSV*"), env(QC_exit), emit: to_minimap2
  tuple val(sampleId), env(TYPE), emit: json
script:
"""
# genome jest w kontenerze

if [ ${QC_STATUS} == "nie" ]; then
  QC_exit="nie"
  TYPE="unk"
  REF_GENOME_ID="unk"
  touch pairs.tsv
  touch RSV_dummy.fasta
  touch RSV_dummy.fasta.amb
  touch primers.bed
else
  REFERENCE_GENOME_FASTA="/home/data/rsv/genome/RSV/RSV.fasta"
  minimap2 -a -x map-ont -t ${params.threads} -o tmp.sam \${REFERENCE_GENOME_FASTA} ${reads}  >> ${log} 2>&1
  samtools view -@ ${params.threads} -Sb -o type_determination.bam -F 2052 tmp.sam

  ILE_A=`samtools view type_determination.bam | grep "hRSV/A/England/397/2017" | wc -l `
  ILE_B=`samtools view type_determination.bam | grep "hRSV/B/England/RE20000104/2020" | wc -l `

  if [[ \${ILE_A} -lt 1000 && \${ILE_B} -lt 1000 ]]; then
      QC_exit="nie"
      TYPE="unk"
      REF_GENOME_ID="unk"
      touch pairs.tsv
      touch RSV_dummy.fasta
      touch RSV_dummy.fasta.amb
      touch primers.bed
  elif [[ \${ILE_A} -gt 1000 || \${ILE_B} -gt 1000 ]]; then
    QC_exit="tak"
    if  [ \${ILE_A} -gt \${ILE_B} ]; then
      TYPE="A"
      cp /home/data/rsv/primers/A/${params.primers_id}/*bed primers.bed
      cp /home/data/rsv/primers/A/${params.primers_id}/*tsv pairs.tsv
      cp /home/data/rsv/genome/RSV_A/RSV* .
    else
      TYPE="B"
      cp /home/data/rsv/primers/B/${params.primers_id}/*bed primers.bed
      cp /home/data/rsv/primers/B/${params.primers_id}/*tsv pairs.tsv
      cp /home/data/rsv/genome/RSV_B/RSV* .

    fi
  REF_GENOME_ID=`head -1 genome.fasta | cut -d " " -f1 | tr -d ">"`
  fi
fi


"""

}
