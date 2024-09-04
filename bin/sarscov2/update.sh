#!/usr/bin/env bash

# This script shouldn't be run directly. It should be run inside the updater container.

# Function to update the nextclade database
update_nextclade() {
    local db_path=$1
    /opt/nextclade/bin/nextclade dataset get --name sars-cov-2 --output-zip "${db_path}/sars-cov-2.zip"
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
    elif [ "$db_name" == "kraken" ]; then
        update_kraken "$db_path"
    elif [ "$db_name" == "freyja" ]; then
        update_freyja "$db_path"
    fi
}

# Check if argument is provided
if [ -z "$1" ]; then
    echo "No argument supplied. Please provide one of the following: nextclade, pangolin, kraken, freyja or all"
    exit 1
fi

# Check if argument is valid
if [ "$1" != "nextclade" ] && [ "$1" != "pangolin" ] && [ "$1" != "kraken" ] && [ "$1" != "freyja" ] && [ "$1" != "all" ]; then
    echo "Invalid argument supplied. Please provide one of the following: nextclade, pangolin, kraken, freyja or all"
    exit 1
fi

# Update all databases if 'all' is specified
if [ "$1" == "all" ]; then
    update_database "nextclade"
    update_database "pangolin"
    update_database "kraken"
    update_database "freyja"
    exit 0
fi

# Update the specified database
update_database "$1"
exit $?