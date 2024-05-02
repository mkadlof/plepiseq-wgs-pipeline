# Pipeline Steps

## Main section
- trimmomatic - removing adapters
- bwa - mapping reads
- filtering - filtering bad quality reads
- masking - masking primers
- merging - intermediate step
- viterbi - improving alignment
- lowCov - detecting low coverage regions
- varScan - small INDEL caller
- freeBayes - small INDEL caller
- lofreq - small INDEL caller
- consensus - consensus sequence of three callers
- consensusMasking - masking low coverage regions

## Quality Check section
- fastqc_1 - check raw reads
- fastqc_2 - check reads after trimmomatic

## Contaminations section
- Kraken2 - detecting other reads from other organisms

## Coinfections section
- coinfections_ivar - intermediate step
- coinfections_varscan - intermediate step
- coinfection_analysis - detecting coinfections 
- freyja - detecting coinfections - alternative method

## Dehumanization section
- dehumanization - removing non viral reads

## Functional analysis
- vcfForFasta - intermediate step
- snpEff - detecting effects of mutations on protein sequence

## SV Calling section
- picard - intermediate step
- manta - detecting large SVs

## Lineages section
- pangolin - detecting pango line
- nextclade - detecting nextclade line

## Model building section
- modeller - building Spike protein model