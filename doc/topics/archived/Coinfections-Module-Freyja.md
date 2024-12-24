# Coinfections: Module Freyja

```Bash
        mkdir variants_files depth_files demix_files
        freyja variants mapped_reads.bam --variants variants_files/test.variants.tsv --depths depth_files/test.depth --ref ${reference_fasta}
        freyja demix variants_files/test.variants.tsv depth_files/test.depth --output demix_files/test.output --confirmedonly --barcodes  /home/external_databases/freyja/usher_barcodes.csv
        freyja aggregate demix_files/ --output coinfections.tsv
```

Freyja is an alternative tool to detect coinfections in SARS-CoV-2 samples.

> The method uses lineage-determining mutational “barcodes” derived from the UShER global phylogenetic tree as a basis set to solve the constrained (unit sum, non-negative) de-mixing problem.

[Documentation](https://andersen-lab.github.io/Freyja/index.html)