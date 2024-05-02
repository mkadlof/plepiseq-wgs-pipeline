# Main: Module lowCov

```Bash
    pysam_quality_mask_final.py ${bam} 10 ${params.mask}
```

To identify regions with low coverage, we utilize the "count coverage" function available within the pysam package. This function has an advantage over a similar function in the `bedtools` package because it allows counting coverage while considering nucleotide quality. The Python script takes 3 arguments. The first is the bam file from which coverage is calculated, the second is the quality threshold, only nucleotides mapping to a region with at least this value are considered. The last argument is the coverage value, only positions in the genome with coverage below this value are returned. The output is a file named `quality_mask.bed`. In this file, each position with coverage below the threshold is returned as a separate entry. The file has standard 3 columns. The first is the chromosome name, the second is the start of the region (0-indexed), and the third column is the end of the region (0-indexed). Like any bed file, the regions are closed on the left and open on the right. An example content of the file is:

```csv
MN908947.3    0    1
MN908947.3    1    2
MN908947.3    2    3
MN908947.3    2    4
MN908947.3    100    200
```

```Bash
    cat quality_mask.bed | bedtools merge -d 2 | \
                           awk 'BEGIN {OFS = "\t"}; {if (\$3-\$2 > 3) print \$1,\$2,\$3}' >> low_coverage.bed
```
In the next step, using `bedtools`, regions in the quality_mask.bed file are merged if the ranges they cover are separated by no more than 2 nucleotides. If, after merging, a low coverage region is at least 4 nucleotides long, it is saved to the `low_coverage.bed file`; otherwise, the region is ignored. This empirical criterion helps avoid situations where coverage in a region oscillates around the value given in the `$mask` argument, and the region is alternately masked and unmasked by very short segments. The example content of the quality_mask.bed file shown above, after this step, would look like this:

```csv
MN908947.3    0    4
```

This information is stored in the `low_coverage.bed` file.

```Bash
    bedtools maskfasta -fi ${reference_fasta} \
                       -bed low_coverage.bed \
                       -fo lowcoverage_masked.fa
```

The final step is to create a fasta file with the reference genome where positions with low coverage are masked with Ns. For this purpose, we use the appropriate function of the `bedtools` program. The fasta file with such a genome is named `lowcoverage_masked.fa`.
