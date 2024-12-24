# Functional analysis: vcfForFasta

```Bash
    prep_own_vcf.py ${reference_fasta} consensus.fa \$N \${vcf_input} \${vcf_output}
```

A simple module designed to transform a fasta file containing a consensus sequence derived from three programs (varScan, FreeBayes, and Lofreq) into a VCF file format. The program has one parameter `-N`, which specifies the maximum gap size between mutations to merge them into a single complex mutation. This parameter is set within the module, and changing it requires rebuilding the container `nf_illumina_sars-3.0-main`.
