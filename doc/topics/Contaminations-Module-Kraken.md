# Contaminations: Module Kraken2

```Bash
    kraken2 --db /home/external_databases/kraken2 \
            --report raport_kraken2.txt \
            --threads ${params.threads} \
            --gzip-compressed \
            --minimum-base-quality 30 \
            --use-names ${reads[0]} ${reads[1]} >> raport_kraken2_individualreads.txt 2>&1

```
Kraken2 is a tool capable of detecting contamination in a sample by aligning reads to genomes in a database and assessing their taxonomic origin. The principle of Kraken involves mapping reads to genomes and performing read-based classification against a reference database.

As a result, users obtain information about the proportions of reads originating from different organisms.

Kraken requires a large database for classification. Details are described in [updates.md#external-databases-updates](updates.md#external-databases-updates).
