# Quickstart

### Installation and usage

1. Install NextFlow

```bash
  curl -s https://get.nextflow.io | bash
  mv nextflow ~/bin
```

2. Build two containers:

```bash
  docker build --target production -f Dockerfile-main -t nf_illumina_sars-3.0-main .
  docker build --target prodcution -f Dockerfile-manta -t nf_illumina_sars-3.0-manta .
```

3. Install `pangolin-data` package locally in `data/pangolin` directory:

```bash
pip install \ 
    --target data/pangolin \
    --upgrade \
    git+https://github.com/cov-lineages/pangolin-data.git@<VERSION>
```

Substitute `<VERSION>` with desired tag name from git repository.
List of tags is available [here](https://github.com/cov-lineages/pangolin-data/tags).

4. Copy `run_nf_pipeline.sh.template` to `run_nf_pipeline.sh` and fill in the paths to the reads and output directory.
5. Run the pipeline:

```bash
  ./run_nf_pipeline.sh
```