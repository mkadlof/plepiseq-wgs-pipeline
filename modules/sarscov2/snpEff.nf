process snpEff {
    tag "snpEff:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}_detected_variants_consensus_annotated.txt"

    input:
    tuple val(sampleId), path(consensus_vcf_gz), path(consensus_vcf_gz_tbi), val(QC_status_vcf), path('forvariants.bam'), path('forvariants.bam.bai'), val(QC_status), path(ref_genome)
    output:
    tuple val(sampleId), path("${sampleId}_detected_variants_consensus_annotated.txt")

    script:
    """
    
    if [[ ${QC_status_vcf} == "nie" && ${QC_status} == "nie" ]]; then
      touch ${sampleId}_detected_variants_consensus_annotated.txt
    else
      
      # Check which organism we are annalyzing ?
      # Why aren't we using params.species ?
      if [ `head -1 ${ref_genome} | awk '{split(\$1, a, "/"); {if (a[2] == "A") {print 1} else {print 0} } }'` == 1 ]; then
        snp_eff='hRSV_A'
        cp /home/data/rsv/snpEff/config/snpEff.config /opt/snpEff/snpEff.config
        cp -r /home/data/rsv/snpEff/data/* /opt/snpEff/data/
      elif  [ `head -1 ${ref_genome} | awk '{split(\$1, a, "/"); {if (a[2] == "B") {print 1} else {print 0} } }'` == 1 ]; then
        snp_eff='hRSV_B' 
        cp /home/data/rsv/snpEff/config/snpEff.config /opt/snpEff/snpEff.config
        cp -r /home/data/rsv/snpEff/data/* /opt/snpEff/data/
      elif  [ `head -1 ${ref_genome} | awk '{if (substr(\$0, 2) == "MN908947.3") {print 1} else {print 0}}'` == 1 ]; then
        snp_eff='MN908947.3'
      else
        snp_eff='GRYPA' # NIE ZAIMPLEMENTOWANE
      fi

      if [ \${snp_eff} == "GRYPA" ] ;then
        ## dummy output for now
        touch ${sampleId}_detected_variants_consensus_annotated.txt
      else 
        # either RSV or SARS, still let us assume these organisms have multiple segments
        zcat ${consensus_vcf_gz} | grep -v "^#" | awk '{print \$1, \$2-1, \$2}' | tr " " "\t" >> tmp.bed
        freebayes -t tmp.bed --fasta-reference ${ref_genome} --min-coverage 1 --min-mapping-quality 1 --min-base-quality 1 --ploidy 1 forvariants.bam >> freebayes_consesus_only.vcf

        java -jar /opt/snpEff/snpEff.jar ann -noStats \${snp_eff} ${consensus_vcf_gz} > detected_variants_consensus_annotated.vcf
        bgzip --force detected_variants_consensus_annotated.vcf
        tabix detected_variants_consensus_annotated.vcf.gz

        # split multi segment reference genome into individual fastas 
        cat ${ref_genome} | awk '{if (substr(\$0, 1, 1)==">") { new_name=\$0; gsub("\\\\.", "_", new_name); gsub("/", "_", new_name);  filename=("reference_"substr(new_name,2) ".fasta"); print \$0 >> filename } else {print toupper(\$0)  >> filename}}'


        for PLIK in `ls reference_*`; do
          SEGMENT=`basename \${PLIK} ".fasta" | cut -d "_" -f2- `
          SEGMENT_EXACT_NAME=`head -1 \${PLIK} | tr -d ">"` # fo extracting data from vcf
        
          # extracting infromations from snpEff output for analyzed segment
          # We nned to include CHROM name in the outpur in this version
          bcftools query --format '%CHROM | %POS | %REF%POS%ALT| %ANN \\n' \
                  detected_variants_consensus_annotated.vcf.gz | \
                  cut -d "|" -f1,2,3,5,7,14 | \
                  tr "|" "\\t" | \
                  awk  'BEGIN {OFS = "\\t"} {if ( \$4 == "upstream_gene_variant" || \$4 == "downstream_gene_variant") {gene="."; aa="."} else {gene=\$5; aa=\$6}; print \$1, \$2, gene, \$3, aa, \$4}'  | grep \${SEGMENT_EXACT_NAME} > part1_\${SEGMENT}.txt
        
          ### extracting DP and allele usage from freebayes output for analyzed segment

          cat freebayes_consesus_only.vcf |  grep "^#\\|\${SEGMENT_EXACT_NAME}" >> freebayes_consesus_only_\${SEGMENT}.vcf
          bcftools query --format "%POS %INFO/DP %INFO/SRF %INFO/SRR %INFO/SAF %INFO/SAR %POS \n" freebayes_consesus_only_\${SEGMENT}.vcf | awk '{print \$7, \$2, (\$6+\$5)/(\$3+\$4+\$5+\$6)}' >> part2_\${SEGMENT}.txt
        
          ### merging two files if freebayes dosn't have info regarding a variant from consensus.vcf we
          ### put "-" in DP and AF columns
          join -a 1 -1 2 -2 1 -o1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.3  -e '-' part1_\${SEGMENT}.txt part2_\${SEGMENT}.txt | sort -nk2 | cut -d " " -f1- | tr " " "\\t" >> detected_variants_consensus_annotated_\${SEGMENT}.txt
        done
        # merge results for all the segments 
        cat detected_variants_consensus_annotated_*.txt >> ${sampleId}_detected_variants_consensus_annotated.txt
      fi # koniec if-a na organizm
    fi # koniec if-a na QC
   """
}
