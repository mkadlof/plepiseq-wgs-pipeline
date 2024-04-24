# External databases updates

Some components of the pipeline require access to their two databases, which are updated roughly once every two weeks. Different software pieces need to be updated in different ways. To make this process as smooth and painless as possible, we prepared a dedicated Docker container exactly for this task, along with a bash script for running it with appropriate volume mounts. The script should be placed in either the cron or systemd timer and run on a weekly basis.

## Updates procedure {id="updates_procedure"}
Build the dedicated container:

```bash
docker build --target updater -f Dockerfile-main -t nf_illumina_sars-3.0-updater:latest .
```

Run the updater script. The working dir must be in project root directory.

```bash
update_external_databases.sh nextclade
update_external_databases.sh pangolin
update_external_databases.sh kraken
update_external_databases.sh freyja
```

Total size of downloads is ~55 GiB.

| Database  | Size     |
|-----------|----------|
| pangolin  | ~90 MiB  |
| nextclade | ~1.3 MiB |
| kraken    | ~55 GiB  |
| freyja    | ~100 MiB |

If everything work fine in directories `data\pangolin` and `data\nextclade` you should see downloaded content like below:

```
data/nextclade/
└── sars-cov-2.zip

data/pangolin/
├── bin
├── pangolin_data
└── pangolin_data-1.25.1.dist-info

data/kraken
└── k2_standard_20240112.tar.gz

data/freyja/
├── curated_lineages.json
├── lineages.yml
└── usher_barcodes.csv
```

It is recommended to put the following in crontab or equivalently systemd timer.

```bash
0 3 * * 6 cd /path/to/sars-illumina && bin/update_external_databases.sh nextclade
5 3 * * 6 cd /path/to/sars-illumina && bin/update_external_databases.sh pangolin
10 3 * * 6 cd /path/to/sars-illumina && bin/update_external_databases.sh freyja
15 3 1 */3 * cd /path/to/sars-illumina && bin/update_external_databases.sh kraken
```

## Updates internals

The following section contain information what and how is updated. Unless you need to debug or refactor the code, and you followed guides in chapter [](updates.md#updates_procedure) you can safely skip it.

### List of components that require updates

* Nextclade
* Pangolin
* Kraken
* Freyja


