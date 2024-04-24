# Updateing Nextclade database

[Nextclade](https://docs.nextstrain.org/projects/nextclade/) is software for assigning evolutionary lineage to SARS-Cov2.
To make it work properly, it requires a database which is updated roughly once every two weeks.

The recommended way of downloading dataset is using `nxtclade` tool.

```bash
nextclade dataset get --name sars-cov-2 --output-zip sars-cov-2.zip
```

Detailed manual is available [](https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html).
Nextclade is downloading `index.json` from site: [](https://data.clades.nextstrain.org/v3/index.json), and based on that files it decide what to download and from where. Probably the same data are available directly on GitHub: [](https://github.com/nextstrain/nextclade_data/tree/master/data/nextstrain/sars-cov-2/wuhan-hu-1/orfs).

The command above will download a `sars-cov-2.zip` file in desired destination (default: `data/nextclade`). That directory have to be mounted inside `main` container. It is done by Nextflow in the `modules/nextclade.nf` module.

```Groovy
process nextclade {
    (...)
    containerOptions "--volume ${params.nextclade_db_absolute_path_on_host}:/home/SARS-CoV2/nextclade_db"
    (...)
```