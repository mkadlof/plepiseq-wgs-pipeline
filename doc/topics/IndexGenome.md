# IndexGenome

```Bash
    samtools faidx ${reference_fasta}
```

Many modules require a reference genome as input (either a standalone FASTA file or with a FAI index). This module creates an index once and passes it to the required modules. For clarity in the pipeline visualization diagram, this module is hidden.
