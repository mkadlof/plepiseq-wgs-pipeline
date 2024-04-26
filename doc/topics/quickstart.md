# Quickstart

### Installation and usage

1. Install [Docker](https://docs.docker.com/desktop/install/linux-install/) 
2. Install NextFlow

```bash
  curl -s https://get.nextflow.io | bash
  mv nextflow ~/bin
```

3. Clone the repository:

```bash
  git clone %git.url%
````  

4. Copy `third-party/modeller/config.py.template` to `third-party/modeller/config.py` and replace the line

```license = 'YOUR_MODELLER_KEY'```

with the actual Modeller key you own. If you don't have one, you can get a free academic license [here](https://salilab.org/modeller/registration.html).

5. Build three containers:

```bash
  docker build --target production -f Dockerfile-main -t nf_illumina_sars-3.0-main .
  docker build --target prodcution -f Dockerfile-manta -t nf_illumina_sars-3.0-manta .
  docker build --target updater -f Dockerfile-main -t nf_illumina_sars-3.0-updater:latest .
```

6. Download latest version of external databases:

In project root dir run:
```bash
  ./update_external_databases.sh
```
This should fill directories in `data/pangolin` and `data/nextclade`.
For more details read the chapter [](updates.md).

7. Copy `run_nf_pipeline.sh.template` to `run_nf_pipeline.sh` and fill in the paths to the reads and output directory.

8. Run the pipeline:

```bash
  ./run_nf_pipeline.sh
```