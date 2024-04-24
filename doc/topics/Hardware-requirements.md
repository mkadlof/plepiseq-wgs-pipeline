# Hardware requirements

## Platform

Pipeline is intended to run on Unix-like computing server with x86_64 architecture.

## Memory and CPU

Pipeline consist of multiple steps (processes) that are run in separate containers possibly concurrently. Each process
has its own hardware requirements. Proceses vary highly in their demands, from very low to very high. Some may benefit
from multiple cores while others are single-threaded or fast enough to not require more than one core. Exact
requirements depend on expected number of samples analyzed in parallel.

During our tests we run the pipeline for 32 samples in parallel on a machine with 96 cores (Intel(R) Xeon(R) Gold 6240R
CPU @ 2.40GHz) and 503 GiB of RAM, and we faced out-of-memory issues. We introduced limits on memory-intensive
processes (BWA and Kraken2) to maximum 5 concurrent instances, to avoid OOM killer.

### Single sample mode

In case of running the pipeline in single sample mode we recommend using at least 16 cores and 64 GiB of RAM.

### Multiple samples mode

In case of running set of samples in parallel we recommended using 4 cores and ~100 GiB of RAM per sample.

## Storage requirements

Total storage requirements is **~60 GiB** of constant data, and further **~3.3 GiB per sample**.

### Performance

Many processes in the perform a lot of I/O operations, thus pipeline definitely can benefit from fast storage. We
recommend store external databases, and temporary files in fast storage like NVMe SSD in RAID 0. It also may be
beneficial to store it on in-memory filesystem like `tmpfs`, however no extensive tests were performed.

### Docker images sizes:

Pipeline consist of three docker images two for computations and one is wrapper for external databases updates.

| Image                          | Size         |
|--------------------------------|--------------|
| `nf_illumina_sars-3.0-main`    | 1.82 GiB     |
| `nf_illumina_sars-3.0-updater` | 257 MiB      |
| `nf_illumina_sars-3.0-manta`   | 1.26 GiB     |
| **Total**                      | **3.33 GiB** | 

### Databases sizes:

Pipeline require access to external databases. Total size of databases is **~56 GiB**.

| Database  | Size        |
|-----------|-------------|
| pangolin  | ~90 MiB     |
| nextclade | ~1.3 MiB    |
| kraken    | ~55 GiB     |
| freyja    | ~100 MiB    |
| **Total** | **~56 GiB** | 

### Temporary files and results

Pipeline generates a lot of temporary files, which are stored in `work` directory.
According to our tests on EQA2023 dataset (32 samples) took ~105 GiB of disk space **~3.3 GiB per sample**.
Please note that those tests are not representative, and real life data may vary significantly.

## GPU requirements

Pipeline does not exploit GPU acceleration, so no GPU is required.