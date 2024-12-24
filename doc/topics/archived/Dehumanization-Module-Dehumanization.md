# Dehumanization: Module Dehumanization

```Bash
    samtools view mapped_reads.bam | cut -f1 | sort | uniq >> lista_id_nohuman.txt
    seqtk subseq ${reads[0]} lista_id_nohuman.txt >> forward_paired_nohuman.fq
    seqtk subseq ${reads[1]} lista_id_nohuman.txt >> reverse_paired_nohuman.fq

    gzip forward_paired_nohuman.fq
    gzip reverse_paired_nohuman.fq
   ```

The goal of this module is to remove reads that do not map to the reference genome to avoid loading human reads into external databases. It creates a new set of purified ang gziped FASTQ files.
