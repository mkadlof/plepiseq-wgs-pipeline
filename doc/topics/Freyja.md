# Updateing Freyja database

[Freyja](https://andersen-lab.github.io/Freyja/index.html) is a tool to recover relative lineage abundances from mixed SARS-CoV-2 samples from a sequencing dataset. 

Natively Freyja updates require installing Freyja python package (by default with conda) and running `freyja update` command. This command will download the latest version of curated lineages and usher barcodes.
However, in our pipeline we do it simpler way. We download the files directly from the dedicated GitHub repository [andersen-lab/Freyja-data](https://github.com/andersen-lab/Freyja-data) using regular `wget`, which is implemented in the `update_external_databases.sh` script, and `updater` container.
