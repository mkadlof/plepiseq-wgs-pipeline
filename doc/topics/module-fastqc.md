# Module: fastqc

## Description
The `fastqc` module performs quality control (QC) on paired-end sequencing reads. It generates statistical summaries and visualizations of read quality, length, and position-based metrics. It also evaluates if the input reads pass predefined QC thresholds, producing JSON files with detailed results and an overall QC status.

## Input
- **sampleId**: Identifier for the sample being analyzed.
- **reads**: Path to paired-end sequencing read files (forward and reverse).
- **QC_STATUS**: Initial QC status from previous module (default: "tak").
- **prefix**: Prefix for output files.

> The **prefix** field is used to identify module calls. The fastqc module is invoked twice: before and after Trimmomatic, and the quality of the reads is assessed twice. The same module is imported twice. To allow Nextflow to distinguish between them, a prefix field is added with an arbitrary label: "pre-filtering" or "post-filtering".
{style="note"}

## Output
- CSV files with read quality metrics and histograms (e.g., read quality histogram, read length histogram, position-based quality plot).
- JSON files (`forward_<prefix>.json`, `reverse_<prefix>.json`) containing detailed QC results.
- Environment variable **QC_STATUS_EXIT** indicating the overall QC status ("tak", "nie", or "blad").

> The QC_STATUS_EXIT field acts as a quality control flag. It is passed to each subsequent module. If an error is detected at any stage that prevents further analysis, or if data quality issues are identified, the flag is set to "error" or "no." In such cases, subsequent modules will not perform analyses but will instead pass the results obtained so far to the result-aggregating module.
{style="note"}