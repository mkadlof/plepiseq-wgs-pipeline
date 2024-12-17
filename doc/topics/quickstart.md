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
  docker build --target main -f Dockerfile-main -t pzh_pipeline_viral-4.1-main .
  docker build --target manta -f Dockerfile-manta -t pzh_pipeline_viral-4.1-manta .
  docker build --target updater -f Dockerfile-main -t nf_illumina_sars-4.1-updater .
```

> **TIP**  
> You may add `--no-cache` to avoid caching effects. 


6. Download latest version of external databases:

In project root dir run:
```bash
  ./update_external_databases.sh pangolin
  ./update_external_databases.sh nextclade
  ./update_external_databases.sh kraken
  ./update_external_databases.sh freyja
```
For more details read the chapter [](updates.md).

7. Copy `run_nf_pipeline.sh.template` to `run_nf_pipeline.sh` and fill in the paths to the reads and output directory.

8. Run the pipeline:

```bash
  ./run_nf_pipeline.sh
```