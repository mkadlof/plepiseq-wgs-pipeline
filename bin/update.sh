#!/usr/bin/env bash

# This script shouldn't be run directly. It should be run inside the updater container.

# Check if argument is provided
if [ -z "$1" ]
then
    echo "No argument supplied. Please provide one of the following: nextclade, pangolin, kraken"
    exit 1
fi

# Check if argument is valid
if [ "$1" != "nextclade" ] && [ "$1" != "pangolin" ] && [ "$1" != "kraken" ]
then
    echo "Invalid argument supplied. Please provide one of the following: nextclade, pangolin, kraken"
    exit 1
fi

# Update nextclade
if [ "$1" == "nextclade" ]
then
  /opt/nextclade/bin/nextclade dataset get --name sars-cov-2 --output-zip /home/nextclade/sars-cov-2.zip
fi

# Update pangolin
if [ "$1" == "pangolin" ]
then
  pip install \
      --target pangolin \
      --upgrade \
      git+https://github.com/cov-lineages/pangolin-data.git
fi

# Update kraken2
if [ "$1" == "kraken" ]
then
  /home/kraken_updater.py standard /home/kraken
fi