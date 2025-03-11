In this directory, we place our own versions of primer schemes. The instructions are described in 
the documentation, point B.3. (check the new documentation!)

The schemes downloaded from the internet have been cleaned by me. That is, the primer name always
starts with "nCoV-2019_". Occasionally, a missing column 6 has been added. The file with the scheme
in the given subdirectory is always named `nCoV-2019.scheme.bed`. Missing `pairs.tsv` files
necessary for `ivar` have been added.

### Directories:

V1 to V5 - schemes downloaded from the repository
[https://github.com/artic-network/primer-schemes/tree/master](https://github.com/artic-network/primer-schemes/tree/master)
standardly used in the artic protocol.

Directories:
- EQA2023.SARS1
- EQA2023.SARS2
- EQA2024.V4_1 (practically identical, the only difference is the range of primer 64_LEFT)
- EQA2024.V5_3 (identical to the directory V.5.3.2 because V.5.3.2 is its copy)

contain primers used in EQA tests according to the names they had there.

VarSkip2 - scheme used by the NEBNext VarSkip Short v2 protocol, downloaded from
[https://github.com/nebiolabs/VarSkip/blob/main/neb_vss2a.primer.bed](https://github.com/nebiolabs/VarSkip/blob/main/neb_vss2a.primer.bed),
has been corrected and adapted to our pipeline. Therefore, it is DIFFERENT than in the repo for
nanopore data.

obserco_extra - directory with primers used in PZH as part of the obserco grant, it is in fact
midnight1200 with added extra primers, but most of these additions have different sequences but map
to the same genomic region as the original primer.

It should be noted that ALL primers in this directory (save EQA2024.V4_1.nanopore) are "artificially" extended so that the
amplicon is 1bp longer in the 5' and 3' directions. This allows `ivar` to properly mask reads that
map just 1bp beyond their amplicon. This additional nucleotide should not biologically occur, but I
attribute this to an Illumina error

Temporarily nanopore primers for sars will have .nanopore extension in directory name, the directory should contain a single .bed file, and the "pairs.tsv". The pairs file is NOT required
required by the pipeline and can be empty, but is required by nexflow modules that are shared between nanopore and illumina

