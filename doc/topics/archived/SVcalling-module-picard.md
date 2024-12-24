# SVs calling: Module picard

```Bash
    java -jar /opt/picard/picard.jar PositionBasedDownsampleSam \
                                    --INPUT ${bam} \
                                    --OUTPUT downsample.bam -F \${NORM}
```

At this stage, we prepare a .bam file that will be used by the Manta program for SV identification. The source .bam file is the one generated in [module bwa](Main-Module-bwa.md), not [module merging](Main-Module-merging.md), because primer masking may result in the identification of additional/longer SVs than actually exist.

Picard was not designed for analyzing data from amplicon-based sequencing. Empirically, consistent results with EQA were obtained when the .bam file contained approximately 200,000 reads. To obtain a .bam file with this number of reads while preserving their genome distribution as in the original file, the `PositionBasedDownsampleSam` function was used. Apart from specifying the input and output file names, the function requires an argument (-F), i.e., a fraction (not a number) of the original number of reads to be retained in the newly created file. In the script, this fraction is symbolically denoted as NORMALIZED($max_number_for_SV) to illustrate that it depends on the value provided by the user for the max_number_for_SV variable, which defaults to 200,000.

It's worth explaining why the `PositionBasedDownsampleSam` function was not used in Step 4. This is because amplicon-based sequencing is characterized by very large coverage gaps between different, often adjacent genome regions. Thus, there are regions with coverage as high as 20,000 next to regions with significantly lower coverage, e.g., 50. In such a situation, it is impossible to reduce the coverage in a region with very high coverage to values like 1,000-5,000 while simultaneously maintaining coverage of 50 for a region with low coverage. After applying the `PositionBasedDownsampleSam` function, the coverage in low coverage regions will be much lower, often below 20, which is an important threshold because it defines regions that are masked as "N". Moreover, it may lead to erroneous SV identification in such regions. Additionally, `PositionBasedDownsampleSam` does not understand that reads come from sequencing a specific amplicon, and some reads are "more important" and carry correct information, while some reads can be filtered out because they are irrelevant.
