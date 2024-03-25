#!/usr/bin/env bash

# This script shouldn't be run directly. It should be run inside the updater container.

# Update nextclade
/opt/nextclade/bin/nextclade dataset get --name sars-cov-2 --output-zip /home/nextclade/sars-cov-2.zip

# Update pangolin
pip install \
    --target pangolin \
    --upgrade \
    git+https://github.com/cov-lineages/pangolin-data.git