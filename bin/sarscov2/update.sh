#!/usr/bin/env bash

# This script shouldn't be run directly. It should be run inside the updater container.

# Function to update the nextclade database
update_nextclade() {
    local db_path=$1
    # SARS
    /opt/nextclade/bin/nextclade dataset get --name='sars-cov-2' --output-zip "${db_path}/sars-cov-2.zip"
    # INFL
    /opt/nextclade/bin/nextclade dataset get --name='flu_h1n1pdm_ha' --output-zip "${db_path}/H1N1_HA.zip"
    /opt/nextclade/bin/nextclade dataset get --name='flu_h1n1pdm_na' --output-zip "${db_path}/H1N1_NA.zip"
    /opt/nextclade/bin/nextclade dataset get --name='flu_h3n2_ha' --output-zip "${db_path}/H3N2_HA.zip"
    /opt/nextclade/bin/nextclade dataset get --name='flu_h3n2_na' --output-zip "${db_path}/H3N2_NA.zip"
    /opt/nextclade/bin/nextclade dataset get --name='flu_vic_ha' --output-zip "${db_path}/Victoria_HA.zip"
    /opt/nextclade/bin/nextclade dataset get --name='flu_vic_na' --output-zip "${db_path}/Victoria_NA.zip"
    /opt/nextclade/bin/nextclade dataset get --name='flu_yam_ha' --output-zip "${db_path}/Yamagata_HA.zip"
    # RSV
    /opt/nextclade/bin/nextclade dataset get --name='rsv_a' --output-zip "${db_path}/RSV_A.zip"
    /opt/nextclade/bin/nextclade dataset get --name='rsv_b' --output-zip "${db_path}/RSV_B.zip"
    return $?
}

# Function to update the pangolin database
update_pangolin() {
    local db_path=$1
    pip install \
        --target "${db_path}" \
        --upgrade \
        git+https://github.com/cov-lineages/pangolin-data.git
    return $?
}

# Function to update the kraken database
update_kraken() {
    local db_path=$1
    /home/kraken_updater.py standard "${db_path}"
    return $?
}

# Function to update the freyja database
update_freyja() {
    local db_path=$1
    cd "${db_path}"
    wget https://raw.githubusercontent.com/andersen-lab/Freyja-data/main/lineages.yml
    wget https://raw.githubusercontent.com/andersen-lab/Freyja-data/main/curated_lineages.json
    wget https://github.com/andersen-lab/Freyja-data/raw/main/usher_barcodes.csv
    return $?
}

# Function to update a specific database
update_database() {
    local db_name=$1
    local db_path="/home/external_databases"

    # Ensure the directory exists
    if [ ! -d "$db_path" ]; then
        echo "Directory $db_path does not exist. Creating it..."
        mkdir -p "$db_path"
    fi

    # Call the appropriate update function
    if [ "$db_name" == "nextclade" ]; then
        update_nextclade "$db_path"
    elif [ "$db_name" == "pangolin" ]; then
        update_pangolin "$db_path"
    elif [ "$db_name" == "kraken2" ]; then
        update_kraken "$db_path"
    elif [ "$db_name" == "freyja" ]; then
        update_freyja "$db_path"
    fi
}

# Check if argument is provided
if [ -z "$1" ]; then
    echo "No argument supplied. Please provide one of the following: nextclade, pangolin, kraken2, freyja or all"
    exit 1
fi

# Check if argument is valid
if [ "$1" != "nextclade" ] && [ "$1" != "pangolin" ] && [ "$1" != "kraken2" ] && [ "$1" != "freyja" ] && [ "$1" != "all" ]; then
    echo "Invalid argument supplied. Please provide one of the following: nextclade, pangolin, kraken2, freyja or all"
    exit 1
fi

# Update all databases if 'all' is specified
if [ "$1" == "all" ]; then
    update_database "nextclade"
    update_database "pangolin"
    update_database "kraken2"
    update_database "freyja"
    exit 0
else
  # Update the specified database
  update_database "$1"
fi

exit $?
