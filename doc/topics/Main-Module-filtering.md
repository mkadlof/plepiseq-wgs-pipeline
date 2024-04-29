# Main: Module filtering

This script is extremely important. It checks a received BAM file and performs quasi-undersampling. The parameters selected above can only be modified by modifying the Python script. The following actions are performed within this script:

1. Identification and saving for further analysis of paired reads located within an amplicon. A paired read is considered to be inside an amplicon if (I) the start of mapping of such a pair is at least at the position corresponding to the first nucleotide of the left primer of the amplicon, and (II) the end of mapping of such a pair is before the last nucleotide of the right primer of the amplicon. This means that for a pair of reads, it is not required for at least one of the reads to map to any of the primers. What is required is not exceeding the boundaries set by a particular amplicon. At the same time, an empirical criterion limiting the number of reads subjected to analysis is applied:
   1. The number of read pairs whose mapping start is between the position of the first nucleotide of the left primer of a given amplicon and the position **up to 150** nucleotides towards the 3’ end cannot exceed **4000**.
   2. The number of reads from point a and reads whose mapping **starts at least 150 nucleotides towards the 3’** end counting from the first nucleotide of the left primer cannot exceed **8000**.
   3. The number of reads from points a and b, as well as reads whose mapping **starts at least 250 nucleotides towards the 3’** end counting from the first nucleotide of the left primer, cannot exceed **12000**.
    Reads that are inside an amplicon but are not selected for further analysis due to exceeding the aforementioned threshold values **are not** further analyzed. Reads that do not meet the criterion of belonging to a single amplicon may undergo further analysis only if low usage is detected for at least one of the amplicons in the analyzed sample. Empirically, usage is considered low if fewer than **100** paired reads are assigned to a given amplicon.
   {type="alpha-lower"}

2. In the second step, reads that most likely belong to a single amplicon, and their source is not a fusion of two neighboring amplicons within the same pool, are removed. Such reads must meet one of two criteria:
    1. They start at least **10** nucleotides towards the **5’ end** counting from the first nucleotide of the left primer, and end before the last nucleotide of the right primer of the amplicon. The amplicon covering such a read must have high usage.
    2. They start at the position corresponding to the first nucleotide of the left primer of the low-usage amplicon and end no further than **10** nucleotides towards the **3’ end** counting from the last nucleotide of the right primer corresponding to the neighboring amplicon. The amplicon covering such a read must have high usage.
    {type="alpha-lower"}

3. In the last step, reads remaining after filtering from points 1 and 3 are analyzed for the possibility of originating from an insert that is a fusion of a low-coverage amplicon and a neighboring amplicon.
    1. If there is an amplicon towards the 5’ end from the low-usage amplicon, we check if the read pair maps to the region defined by the left primer of the low-usage amplicon towards the 5’ end and the right primer of the neighboring amplicon. If in such a pair 50% of the length of any of the read pairs maps to the low-usage amplicon, and 10% of the length of the other read pair maps to the low-usage amplicon, such a read is included in the analysis. If nucleotides from such a pair cover the left primer of the low-usage amplicon towards the 5’ end from the low-usage amplicon and the right primer of the low-usage amplicon, such positions are MASKED using the ivar program. For this purpose, an ad hoc .bed file is created.
    2. If there is an amplicon towards the 3’ end from the low-usage amplicon, we check if the read pair maps to the region defined by the left primer of the low-usage amplicon and the right primer of the neighboring amplicon towards the 5’ end. If in such a pair 50% of the length of any of the read pairs maps to the low-usage amplicon, and 10% of the length of the other read pair maps to the low-usage amplicon, such a read is included in the analysis. If nucleotides from such a pair cover the left primer of the low-usage amplicon and the right primer of the neighboring amplicon towards the 5’ end, such positions are MASKED using the ivar program. For this purpose, an ad hoc .bed file is created.
   {type="alpha-lower"}
Below are illustrated several scenarios. It is assumed that amplicon 1 has high usage (over 100 mapping read pairs), and amplicon 2 has low usage (below 100). Amplicon 1 is located between P1_L and P1_R; amplicon 2 is located between P2_L and P2_R. Both amplicons come from the same pool. In the schema, "X" denotes any nucleotide. The directionality of reads in the pair is not shown.

Positions of elements in the genome.
```
----P1_L----(amplicon1)----P1_R----P2_L------(amplicon2)------P2_R--
```

Hypothetical read pairs:

(Version 1)
```
---------------XXXXXXXXXXXXXXXXXXXXXXXXX--------------------------
------------------------------XXXXXXXXXXXXXXXXXXXXXXXX------------
```
(Version 2)
```
----------XXXXXXXXXXXXXX------------------------------------------
--------------------------------XXXXXXXXXXXXX---------------------
```
(Version 3)
```
----------XXXXXXXXXXXXXXXXX--------------------------------------
-----------------------------------------------XXXXXXXXXXXXXXXXX-
```
(Version 4)
```
XXXXXXXXXXX-------------------------------------------------------
----------------------XXXXXXXXXXXX--------------------------------
```

Version 1.

The read pair is inside an insert formed by the fusion of amplicons 1 and 2. At least half of the nucleotides from read 2 of the pair cover amplicon 2, and at least 10% of the nucleotides from read 1 of the pair cover amplicon 2. Such a read is included in further analysis. The reads do not cover either P1_L or P2_R, so the read passes without masking.

Version 2.

The read pair is inside an insert formed by the fusion of amplicons 1 and 2. Less than 50% of the nucleotides from both read 1 and read 2 of the pair cover amplicon 2. Such a read is not analyzed further.

Version 3.

The read pair is inside an insert formed by amplicons 1 and 2. At least half of the nucleotides from read 2 of the pair cover amplicon 2, and at least 10% of the nucleotides from read 1 of the pair cover amplicon 2. Such a read is included in further analysis. Read 2 covers P2_R, so the 3’ end of this read (bolded) will be masked.

Version 4.

The read pair is not inside an insert formed by amplicons 1 and 2. The read is not analyzed.