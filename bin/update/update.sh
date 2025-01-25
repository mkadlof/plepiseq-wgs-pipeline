#!/usr/bin/env bash

# Script is intended to run with a docker container 
# hence all paths are HARDCODED
# all subsripts intended to download/update a specific database are located in /home/update
# the "top-level" directory where we place subdirectories ar in /home/external_databases


# If updating/downloading database requires 2-3 lines of code the function is inside this file
# more complex downloads/updates are kept in separate files in /home/updates

# We TRY to distinguish updating vs. downloading functionalities
# This is usually based on the presence or absence of specific files in a given database subdirectory

# This script is intended as a updater for both "viral" and "bacterial" pipelines


# This script understands 3 positional arguments (database name, type of kraken database, bacterial genus)
# All values passed to this script are evaluated by update_external_databases.sh

#############################################
# Function to update the nextclade database #
#############################################


## Nextclade
### No differenc in Updating/Downloading  if /home/external_databases/nextclade directory exist, everything inside it will be removed

update_nextclade() {
    if [ -d "/home/external_databases/nextclade" ]; then
	rm -rf /home/external_databases/nextclade/*
    else
        mkdir /home/external_databases/nextclade
    fi

    db_path="/home/external_databases/nextclade"
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
    cd 
    
    return $?
}

## Pangolin
## No differenc in Updating/Downloading  if /home/external_databases/pangolin directory exist, everything inside it will be removed
## --upgrade option in pip had no effect
update_pangolin() {
    if [ -d "/home/external_databases/pangolin" ]; then
        rm -rf /home/external_databases/pangolin/*
    else 
        mkdir /home/external_databases/pangolin
    
    fi

    pip install --target "/home/external_databases/pangolin" git+https://github.com/cov-lineages/pangolin-data.git
    return $?
}

## Kraken2
## For kraken2 there has an update procedure so we do not remove all the files from /home/external_databases/kraken2

update_kraken2() {
    local kraken2_type=$1
    if [ ! -d "/home/external_databases/kraken2" ]; then
	    mkdir /home/external_databases/kraken2
    fi

    # The script accepts two positional arguments path where database will be placed and type of kraken2 database
    python /home/update/kraken_updater.py /home/external_databases/kraken2 "${kraken2_type}"

    return $?
}

# Freyja
## No update procedure so we remove data if directory is present 
update_freyja() {
    if [ -d "/home/external_databases/freyja" ]; then
	    rm -rf /home/external_databases/freyja/*
    else
	    mkdir /home/external_databases/freyja
    fi

    cd /home/external_databases/freyja

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

#############
# Main code *
#############


db_name=$1
kraken_type=$2
genus=$3
if [ ${db_name} == "all" ];then
        update_kraken2 "$kraken_type" >> /dev/null 2>&1
        update_pangolin >> /dev/null 2>&1
        update_freyja >> /dev/null 2>&1
        update_nextclade >> /dev/null 2>&1
elif [ ${db_name} == "kraken2" ]; then
	update_kraken2 "$kraken_type" >> /dev/null 2>&1
elif [ ${db_name} == "pangolin" ]; then
	update_pangolin >> /dev/null 2>&1
elif [ ${db_name} == "freyja" ]; then
	update_freyja >> /dev/null 2>&1
elif [ ${db_name} == "nextclade" ]; then
	update_nextclade >> /dev/null 2>&1
fi
