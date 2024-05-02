# Pipeline customization

All files required by the module for SNP/INDEL identification are located in the directory `data/generic` if the installation was conducted exactly as described in [](quickstart.md). During the image building process, its contents are copied into the image. To use custom files, you need to appropriately modify the contents of the subdirectories "adapters", "contaminations", "genome", "modeller", "primers" and "vcf_template" before building the image.

## Adapters file

The adapter sequences are located in the subdirectory `adapters`. In cases where the length of the insert (the DNA fragment to be sequenced) is shorter than the number of cycles in sequencing, it may happen that residues of adapter sequences are present at the 3' or 5' end next to the actual insert sequence. This unnecessary part of the read should be trimmed, which requires specifying the adapter sequences used during sequencing. Besides the adapters currently distributed with the trimmomatic program, you can create your own adapter sequence files provided they are in fasta format. Path to `adapters` dir is passed as a one of a parameter. The pipeline assumes the use of adapters placed in the file `TruSeq3-PE-2.fa`.

## Indexed genome

The indexed genome of SARS-CoV-2 is located in the subdirectory `data/generic/genome/SarsCov2`. After building the image, files from this directory are available inside the container in the directory /SARS-CoV2/genome/SarsCov2/. To propose your own version of the genome, follow these steps:

1. Download the genome sequences in fasta format (for example, for the SARS-CoV-2 virus, the reference sequence is available on the website https://www.ncbi.nlm.nih.gov/sars-cov-2/). Save the downloaded file with the name sarscov2.fasta in any directory.
2. Index the genome using the program BWA. Navigate to the directory where you saved the sarscov2.fasta file and execute the following command. Assuming the `bwa` program is in the `$PATH`.

       bwa index sarscov2.fasta

3. Index the genome using the program faidx. Navigate to the directory where you saved the sarscov2.fasta file and execute the following command. Assuming the `samtools` program is in the `$PATH`.

       samtools faidx sarscov2.fasta

4. Copy the contents of the directory with your own indexed genome to `data/generic/genome/SarsCov2`. This way, you overwrite the existing files there.

5. Build the image. The built image will only contain the new genome. The originally used version of the genome will not be available to the container.

6. Note that functional mutation effect predictions will only work if the genome sequence has the header `MN908947.3`. To change this, modify the Dockerfile in the section starting with "#SnpEFF". Add at the end of this section, after the line `RUN java -jar snpEff.jar download MN908947.3`, your own identical command but replacing the name `MN908947.3` with the header used by your new genome. The snpEFF database is regularly updated, but it is possible that our sequence is not included in this database.

7. Your own genome also requires creating your own file with the location of primers if the genome has a different header than `MN908947.3`. Changing the genome sequence may require updating the location of primer sequences in the genome.

## Primers

1. For the SARS-CoV-2 virus, two protocols based on amplicons are currently used: the dominant ARTIC protocol and the long-read Midnight protocol. Over time, new versions and modifications related to emerging virus variants appear. The primers available in May 2022 are located in the "primers" subdirectory and are accessible in the container at the path `data/generic/primers`. Each of these directories contains subdirectories `V1`, `V2`, `V3`, `V4`, `V4.1` (primers used in subsequent versions of the ARTIC protocol), and `V1200` (primers used in the Midnight protocol). Additionally, primers used in the EQA test from April 2023 are added (`SARS1_partmerge_exp` for sequencing from the "SARS1" test and `SARS2_partmerge_exp` for sequencing samples in the "SARS2" test). In each of these directories, there are two files: `nCoV-2019.scheme.bed` and `pairs.tsv`. `nCoV-2019.scheme.bed` is a file containing information about the primer locations in the reference genome, and pairs.tsv contains information about which primer pairs flank each of the resulting amplicons. **When invoking the pipeline, the path to the selected bed file with primers must be provided. There is no "default" version.**

2. When creating custom .bed files, remember that (I) the primer positions in the bed file are indexed from 0, not from 1. (II) The ranges are half-open, meaning the start position is treated as the first position in the genome that contains the primer, and the END position is treated as the first position in the genome that the primer does not cover. Below is an example of how a properly formatted bed file should look:

    ```Plain Text
    MN908947.3 29 54 nCoV-2019_1_LEFT 1 +
    MN908947.3 385 411 nCoV-2019_1_RIGHT 1 -
    MN908947.3 319 342 nCoV-2019_2_LEFT 2 +
    MN908947.3 322 333 nCoV-2019_2_LEFT_alt 2 +
    ```

3. The .bed file must contain the following columns separated by tabs: 
   1. Reference genome name. In case of using a genome other than MN908947.3 (point B.2), always create a primer scheme from scratch.
   2. Start position for the primer
   3. End position for the primer
   4. Primer name. The name is built according to the scheme in which consecutive elements are separated by "_": nCoV-2019_(amplikon number)_(from which side of the amplikon the primer is located, possible expressions are LEFT or RIGHT)_(optional field where we provide the expression 'alt' if a given amplikon has more than one primer flanking it from the same side). For example, `nCoV-2019_2_LEFT_alt` means it is the second primer flanking from the 5' side of amplikon number 2.
   5. Pool identifier from which the amplikon originates. It can appear as a single digit, e.g., "1" or "2", or as a unique text, e.g., `nCoV-2019_1` or `nCoV-2019_2`.
   6. The direction of the strand to which the primer hybridizes.
 
   Publicly available files often have column 7 with the primer sequence, but it is not required. "Basic" and "alt" primers are NOT combined unless we understand the implications for further analysis. Primers can be combined if they are duplicates, meaning they have identical values in columns 2 and 3.

4. The pairs.tsv file should contain only two columns separated by tabs:
   i. Name of the primer flanking the amplicon from the 5' side. The name is identical to column 4 of the .bed file.
   ii. Name of the primer flanking the amplicon from the 3' side. The name is identical to column 4 of the .bed file. 
   If there are more than one primer flanking the amplicon from the same side, provide the primer that generates the longer amplicon.

5. After creating both files, place them in a directory with a unique name inside `data/generic/primers`.

6. Primers used in the Midnight protocol (referred to as `V1200`) are available as an archive on the website [](https://zenodo.org/record/3897530#.Xv5EFpMzadY). These files were prepared based on the protocol described at [](https://www.protocols.io/view/sars-cov2-genome-sequencing-protocol-1200bp-amplic-rm7vz8q64vx1/v6?step=20). Download these files by executing the following commands:

       cd ${HOME}/my_primers/nCoV-2019
       wget https://zenodo.org/record/3897530/files/1200bp_amplicon_bed.tar.gz
       tar -zxf 1200bp_amplicon_bed.tar.gz

7. Primers used in the ARTIC scheme are available in the repository [](https://github.com/artic-network/fieldbioinformatics.git)

8. Due to the operation of the primer masking program, it is recommended to extend the amplicon ranges by 1bp in the 3' and 5' directions. Practically, this means that in the primer file, the START field of the LEFT primer will be one less than implied by the sequence, and the END field of the RIGHT primer will be one more.

