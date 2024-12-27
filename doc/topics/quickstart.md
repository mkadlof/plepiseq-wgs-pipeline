# Quickstart

### Installation and usage

1. Install [Docker](https://docs.docker.com/desktop/install/linux-install/) 
2. Install NextFlow

```bash
  curl -s https://get.nextflow.io | bash
```

3. Put the Nextflow binary in the PATH. E.g.:

```bash
sudo mv nextflow /usr/local/bin
```

4. Clone the repository:

```bash
  git clone --depth 1 %git.url%
````  

5. Enter downloaded repo:

```bash
cd %git.name%
```

6. Copy `third-party/modeller/config.py.template` to `third-party/modeller/config.py` and replace the line

```license = 'YOUR_MODELLER_KEY'```

with the actual Modeller key you own. If you don't have one, you can get a free academic license [here](https://salilab.org/modeller/registration.html).

7. Build three containers:

```bash
  docker build --target main -f Dockerfile-main -t pzh_pipeline_viral-4.1-main .
  docker build --target manta -f Dockerfile-manta -t pzh_pipeline_viral-4.1-manta .
  docker build --target updater -f Dockerfile-main -t nf_illumina_sars-4.1-updater .
```
 
> You may add `--no-cache` flag to avoid caching effects. 

> If you encounter a `certificate verify failed` error during the build process, it may be due to being on a corporate network that injects its own certificate. In this case, add the following flag to the build command, where you can pass the certificate provided by your administrator into the container.
> `--build-arg CERT_FILE="$(cat corporate-certificate.crt)"`

8. Download latest version of external databases:

In project root dir run:
```bash
  ./update_external_databases.sh --database pangolin --output-path /path/to/external_databases 
  ./update_external_databases.sh --database nextclade --output-path /path/to/external_databases 
  ./update_external_databases.sh --database kraken2 --output-path /path/to/external_databases --kraken-type standard_08gb
  ./update_external_databases.sh --database freyja --output-path /path/to/external_databases 
```

> kraken2 require providing database name:
> One of following: 
> "standard" "standard_08gb" "standard_16gb" "viral" "minusb" "pluspf" "pluspf_08gb" "pluspf_16gb" "pluspfp" "pluspfp_08gb" "pluspfp_16gb" "nt" "eupathdb48"

For more details read the chapter [](updates.md).

9. Copy `run_nf_pipeline.sh.template` to `run_nf_pipeline.sh` and fill in the paths to the reads and output directory.

10. Run the pipeline:

```bash
  ./run_nf_pipeline.sh
```