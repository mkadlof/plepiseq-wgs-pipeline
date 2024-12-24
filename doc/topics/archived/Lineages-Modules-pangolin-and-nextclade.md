# Lineages: Modules pangolin and nextclade

```Bash
    pangolin --outfile pangolin_lineage.csv \
             --threads ${params.threads} \
             ${consensus_masked_sv_fa}
             
   nextclade run --input-dataset /home/external_databases/nextclade_db/sars-cov-2.zip \
              --output-csv nextstrain_lineage.csv \
              --output-all nextclade_lineages \
              ${consensus_masked_sv_fa}          
```

After obtaining the genome sequence of the sampled individual, each of the variants (individual sequences from each program as well as the consensus) is analyzed for classification into specific variants. For this purpose, we use the pango and nextclade programs. Their invocation does not require special arguments, as their configuration is done at the container creation level.
