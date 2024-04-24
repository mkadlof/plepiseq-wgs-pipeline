# Updateing Kraken database

[Kraken 2](https://ccb.jhu.edu/software/kraken2/) is a taxonomic classification system using exact k-mer matches.
The pipeline utilizes it to detect contamination in samples. It requires a database of approximately 55 GiB, which could potentially be reduced to 16 GiB or 8 GiB, albeit at the cost of sensitivity and accuracy. It updated roughly quarterly. Skipping updates may result in skipping newer taxa.

Database may be built for user own (not recommended) or be downloaded from aws s3. DB is maintained by Kraken 2 maintainers. Here you can find more [here](https://github.com/BenLangmead/aws-indexes), and [here](https://benlangmead.github.io/aws-indexes/).

Downloading is via aws cli (`apt install awscli`), python library `boto3` or HTTP protocol.
There are several types of databases, that differ with set of organism. We use and recommend using `standard` db, which is quite complete. There is also `nt` db which contain all RefSeq and GenBank sequences, but it's size and processing time is too high for routine surveillance.

Links for http downloads are available [here](https://benlangmead.github.io/aws-indexes/k2).
They are in form of:
`https://genome-idx.s3.amazonaws.com/kraken/k2_DBNAME_YYYMMDD.tar.gz`

where DB name is one of following:
`standard`, `standard_08gb`, `standard_16gb`, `viral`, `minusb`, `pluspf`, `pluspf_08gb`, `pluspf_16gb`, `pluspfp`, `pluspfp_08gb`, `pluspfp_16gb`, `nt`, `eupathdb48`.

Refer to official docs for more details.

In case of our pipeline we use python script that is using `boto3` library. It is part of `nf_illumina_sars-3.0-updater` container. For running this script simply pass `kraken` to the container during running.

We recommend using for this dedicated script: `update_externeal_databases.sh`.

Kraken DB will not be updated if local path already contain the file with the same name.
