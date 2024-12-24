# Modules: kraken2 (illumina and nanopore versions)

## kraken2_illumina

### Description
The `kraken2_illumina` module combines Kraken2 with a species identification process for Illumina reads. It evaluates contamination and updates the QC status based on the presence of reads from the expected genus.

### Input
- **sampleId**: Sample identifier.
- **reads**: Path to paired-end sequencing read files.
- **QC_STATUS**: QC status from the upstream module.
- **EXPECTED_GENUS**: Expected genus for the sample.

### Output
- JSON file (`contaminations.json`) with contamination results.
- **QC_status_contaminations**: Updated QC status as an environment variable.
- **FINAL_GENUS**: Most likely genus as an environment variable.

---

## kraken2_nanopore

### Description
The `kraken2_nanopore` module is similar to `kraken2_illumina` but processes single-end reads from Nanopore sequencing.

### Input
- **sampleId**: Sample identifier.
- **reads**: Path to single-end sequencing read file.
- **QC_STATUS**: QC status from the upstream module.
- **EXPECTED_GENUS**: Expected genus for the sample.

### Output
- JSON file (`contaminations.json`) with contamination results.
- **QC_status_contaminations**: Updated QC status as an environment variable.
- **FINAL_GENUS**: Most likely genus as an environment variable.
