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
    local kraken_type=$2
    echo /home/kraken_updater.py "${db_path}" "$kraken_type"
    /home/kraken_updater.py "${db_path}" "$kraken_type"
    return $?
}

# Function to update the freyja database
update_freyja() {
    local db_path=$1
    cd "${db_path}"
    mkdir H1N1 H3N2 H5Nx FLU-B-VIC RSVa RSVb sarscov2
    
    cd sarscov2
      wget https://raw.githubusercontent.com/andersen-lab/Freyja-data/main/lineages.yml
      wget https://raw.githubusercontent.com/andersen-lab/Freyja-data/main/curated_lineages.json
      wget https://github.com/andersen-lab/Freyja-data/raw/main/usher_barcodes.csv
    cd ..

    ALL_SPECIES=(H1N1 H3N2 H5Nx FLU-B-VIC RSVa RSVb)
    for SPECIES in ${ALL_SPECIES[@]}
    do
	cd ${SPECIES}
          wget https://raw.githubusercontent.com/andersen-lab/Freyja-barcodes/refs/heads/main/${SPECIES}/latest/barcode.csv
          wget https://raw.githubusercontent.com/andersen-lab/Freyja-barcodes/refs/heads/main/${SPECIES}/latest/reference.fasta
          wget https://raw.githubusercontent.com/andersen-lab/Freyja-barcodes/refs/heads/main/${SPECIES}/latest/auspice_tree.json
	cd ..
    done 

    mv FLU-B-VIC Victoria
    mv RSVa RSV_A
    mv RSVb RSV_B
    # Other organisms available in Freyja, but not required in PlEpiSeq.
    # wget -O MEASLESN450_barcodes.csv https://raw.githubusercontent.com/andersen-lab/Freyja-barcodes/refs/heads/main/MEASLESN450/latest/barcode.csv
    # wget -O mpox_barcodes.csv https://raw.githubusercontent.com/andersen-lab/Freyja-barcodes/refs/heads/main/MPX/latest/barcode.csv
    # wget -O MEASLESwholegenome_barcodes.csv https://raw.githubusercontent.com/andersen-lab/Freyja-barcodes/refs/heads/main/MEASLESgenome/latest/barcode.csv
    return $?
}

# Function to update a specific database
update_database() {
    local db_name=$1
    local kraken_type=$2
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
        update_kraken "$db_path" "$kraken_type"
    elif [ "$db_name" == "freyja" ]; then
        update_freyja "$db_path"
    fi
}

# Check if argument is provided
if [ -z "$1" ]; then
    echo "No argument supplied. Please provide one of the following: nextclade, pangolin, kraken2, freyja"
    exit 1
fi

# Check if argument is valid
if [ "$1" != "nextclade" ] && [ "$1" != "pangolin" ] && [ "$1" != "kraken2" ] && [ "$1" != "kraken2" ] && [ "$1" != "freyja" ]; then
    echo "Invalid argument supplied. Please provide one of the following: nextclade, pangolin, kraken2, freyja"
    exit 1
fi

if [ "$1" == "kraken2" ]; then

    VALID_KRAKEN2_DB=("standard" "standard_08gb" "standard_16gb" "viral" "minusb" "pluspf" "pluspf_08gb" "pluspf_16gb" "pluspfp" "pluspfp_08gb" "pluspfp_16gb" "nt" "eupathdb48")
    # Check if database name is provided for kraken2
    if [ -z "$2" ]; then
        echo "No database name supplied for kraken2. Please provide one of the following: ${VALID_KRAKEN2_DB[*]}"
        exit 1
    fi

    # Validate the provided database name
    if [[ ! " ${VALID_KRAKEN2_DB[*]} " =~ " $2 " ]]; then
        echo "Invalid database name for kraken2. Please provide one of the following: ${VALID_KRAKEN2_DB[*]}"
        exit 1
    fi
fi


if [  "$1" == "kraken2" ]; then
    echo update_database "$1" "$2"
    update_database "$1" "$2"
else
    # Update the specified database
    update_database "$1"
fi

exit $?
