# Hardware requirements

## Platform

Pipeline is intended to run on any general purpose GNU/Linux computing server with x86_64
architecture.

## CPU

The pipeline consists of multiple steps (processes) that run in separate containers, potentially
concurrently. Each process has its own hardware requirements, which can vary significantly—from
very low to very high demands. Some processes benefit from multiple cores, while others are
single-threaded or fast enough to not require more than one core. The exact requirements depend on
the number of samples expected to be analyzed in parallel. In general, the pipeline can run with
any number of cores, and additional cores may reduce computation time, especially when analyzing
multiple samples simultaneously.

### Single sample mode

In case of running the pipeline in single sample mode we recommend using at least **8 cores per
sample**.

### Multiple samples mode

In case of running set of samples in parallel we recommended using **4 cores per sample**.

## Memory

The most memory-intensive process is Kraken2, which must load its entire database into memory. The
size of the database determines the overall memory requirement. By default, we use the standard
database, which is approximately 80 GiB. Additional memory is required to support other processes
and the operating system, resulting in a total memory requirement of at least 82 GiB per sample.

> Memory requirements can be significantly reduced by choosing a smaller database (e.g., one
> containing only specific branches of the Tree of Life).
>
> A list of available databases, their contents, and sizes is provided here
> [](https://benlangmead.github.io/aws-indexes/k2).
> {style="note"}

## Storage requirements

Total storage requirements is **~60 GiB** of constant data, and further **~3.3 GiB per sample**.

### Performance

Many processes in the perform a lot of I/O operations, thus pipeline definitely can benefit from
fast storage. We
recommend store external databases, and temporary files in fast storage like NVMe SSD in RAID 0. It
also may be
beneficial to store it on in-memory filesystem like `tmpfs`, however no extensive tests were
performed.

### Docker images sizes:

Pipeline consist of three docker images two for computations and one is wrapper for external
databases updates.

| Image                          | Size         |
|--------------------------------|--------------|
| `nf_illumina_sars-3.0-main`    | 2.01 GiB     |
| `nf_illumina_sars-3.0-updater` | 258 MiB      |
| `nf_illumina_sars-3.0-manta`   | 1.72 GiB     |
| **Total**                      | **3.98 GiB** | 

### Databases sizes:

Pipeline require access to external databases. Total size of databases is **%db.total.size%**.

| Database  | Size                |
|-----------|---------------------|
| pangolin  | %db.pangolin.size%  |
| nextclade | %db.nextclade.size% |
| kraken    | %db.kraken.size%    |
| freyja    | %db.freyja.size%    |
| **Total** | **%db.total.size%** | 

> Storage requirements of databases can be significantly reduced the same way as RAM memory
> requirement.
> {style="note"}

### Temporary files and results

The pipeline generates a large number of temporary files, which are stored in the work directory,
and the actual results data in the results directory. The majority of the space in the results
directory is occupied by dehumanized FASTA files.

The results of our tests are summarized in the table below; however, please note that these data
are not representative of real-life sequencing scenarios, and actual resource demands may vary
significantly.

| Pipeline      | Temporary Storage <br/>Average size per sample. | Results<br/>Average size per sample. |
|---------------|-------------------------------------------------|--------------------------------------|
| ILLUMINA_SARS | 4.43 GiB                                        | 355.41 MiB                           |
| ILLUMINA_INFL | 2.28 GiB                                        | ???.??                               |
| ILLUMINA_RSV  | 1.35 GiB                                        | 117.5 MiB                            |
| NANOPORE_SARS |                                                 |                                      |
| NANOPORE_INFL |                                                 |                                      |
| NANOPORE_RSV  |                                                 |                                      |

Average results per sample.

## GPU requirements

Pipeline does not exploit GPU acceleration, so no GPU is required.