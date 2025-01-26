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


# AMRfinder_plus
## Database has an update mechanism
update_amrfinder() {
	if [ ! -d "/home/external_databases/amrfinder_plus" ]; then
		 mkdir /home/external_databases/amrfinder_plus
	fi
	
	/home/update/download_amrfinder.sh 

}

# Metaphlan
## Database has an update mechanism and dedicated script
update_metaphlan() {
        if [ ! -d "/home/external_databases/metaphlan" ]; then
                 mkdir /home/external_databases/metaphlan
        fi

        /home/update/download_metaphlan.sh

}

# Kmerfinder
## There is no update mechanism so we remove all the files 
update_kmerfinder() {
if [ -d "/home/external_databases/kmerfinder" ]; then
        rm -rf /home/external_databases/kmerfinder/*
else
	mkdir /home/external_databases/kmerfinder/
fi
cd /tmp
git clone https://bitbucket.org/genomicepidemiology/kmerfinder_db.git
cd kmerfinder_db/
bash INSTALL.sh /home/external_databases/kmerfinder/ bacteria

}

# Generic CGE updater
## No update mechanism for CGE 
update_cge() {
	local db=$1
	if [ -d "/home/external_databases/${db}" ]; then
		rm -rf /home/external_databases/${db}/*
	else
		mkdir /home/external_databases/${db}/
	fi
	cd /home/external_databases
	git clone https://bitbucket.org/genomicepidemiology/${db}/
	cd ${db}
	python3 INSTALL.py /home/kma/kma_index non_interactive
}

#VFDB
## No update
### For now the python script has HARDCODED usage of 96 CPUs !!!
update_vfdb() {
        if [ -d "/home/external_databases/vfdb" ]; then
                rm -rf /home/external_databases/vfdb/*
        else
                mkdir /home/external_databases/vfdb/
        fi
        cd /home/external_databases/vfdb
        python3 /home/update/download_vfdb.py
}

# MLST data
## No update mechanism
## Data for different genuses originate from eirther pubmlst or enterobase

update_mlst_campylo() {
# subfunction for update_mlst
	local directory=${1}
	local spec=${2}
	local db=${3}
	if [ -d ${directory} ]; then
		rm ${directory}/*
	else
		mkdir -p ${directory}
	fi
	
	cd ${directory}
	python3 /home/update/download_mlst_campylobacter.py "${spec}" "${db}"
}

update_mlst() {
	local genus=${1}
	# Salmonella Escherichia Campylobacter
	if [ ${genus} == "Campylobacter" ]; then
		if [ -d "/home/external_databases/mlst/Campylobacter" ]; then
			rm -rf /home/external_databases/mlst/Campylobacter/*
		else
			mkdir -p /home/external_databases/mlst/Campylobacter
		fi
		#each species has its own seprate MLST
		SPEC="pubmlst_campylobacter_nonjejuni_seqdef"
		# fetus
		update_mlst_campylo "/home/external_databases/mlst/Campylobacter/fetus" "${SPEC}" "C. fetus MLST"
		# helveticus
		update_mlst_campylo "/home/external_databases/mlst/Campylobacter/helveticus" "${SPEC}" "C. helveticus MLST"
		# concisus
		update_mlst_campylo "/home/external_databases/mlst/Campylobacter/concisus" "${SPEC}" "C. concisus/curvus MLST"
		# hyointestinalis
		update_mlst_campylo "/home/external_databases/mlst/Campylobacter/hyointestinalis" "${SPEC}" "C. hyointestinalis MLST"
		# upsaliensis
		update_mlst_campylo "/home/external_databases/mlst/Campylobacter/upsaliensis" "${SPEC}" "C. upsaliensis MLST"
		# lari
		update_mlst_campylo "/home/external_databases/mlst/Campylobacter/lari" "${SPEC}" "C. lari MLST"
		# insulaenigrae
		update_mlst_campylo "/home/external_databases/mlst/Campylobacter/insulaenigrae" "${SPEC}" "C. insulaenigrae MLST"
		# lanienae
		update_mlst_campylo "/home/external_databases/mlst/Campylobacter/lanienae" "${SPEC}" "C. lanienae MLST"
		# sputorum
		update_mlst_campylo "/home/external_databases/mlst/Campylobacter/sputorum" "${SPEC}" "C. sputorum MLST"
		SPEC="pubmlst_campylobacter_seqdef"
		# jejuni
		update_mlst_campylo "/home/external_databases/mlst/Campylobacter/jejuni" "${SPEC}" "MLST"
	elif [ ${genus} == "Salmonella" ]; then
		if [ -d "/home/external_databases/mlst/Salmonella" ]; then
                        rm -rf /home/external_databases/mlst/Salmonella/*
                else
                        mkdir -p /home/external_databases/mlst/Salmonella
                fi
		cd /home/external_databases/mlst/Salmonella
		python3 /home/update/download_salmonella_mlst.py
	elif [ ${genus} == "Escherichia" ]; then
		if [ -d "/home/external_databases/mlst/Escherichia" ]; then
                        rm -rf /home/external_databases/mlst/Escherichia/*
                else
                        mkdir -p /home/external_databases/mlst/Escherichia
                fi
                cd /home/external_databases/mlst/Escherichia
                python3 /home/update/download_escherichia_mlst.py
	elif [ ${genus} == "all" ]; then
                if [ -d "/home/external_databases/mlst/Salmonella" ]; then
                        rm -rf /home/external_databases/mlst/Salmonella/*
                else
                        mkdir -p /home/external_databases/mlst/Salmonella
                fi
                cd /home/external_databases/mlst/Salmonella
                python3 /home/update/download_salmonella_mlst.py

		if [ -d "/home/external_databases/mlst/Escherichia" ]; then
                        rm -rf /home/external_databases/mlst/Escherichia/*
                else
                        mkdir -p /home/external_databases/mlst/Escherichia
                fi
                cd /home/external_databases/mlst/Escherichia
                python3 /home/update/download_escherichia_mlst.py

		if [ -d "/home/external_databases/mlst/Campylobacter" ]; then
                        rm -rf /home/external_databases/mlst/Campylobacter/*
                else
                        mkdir -p /home/external_databases/mlst/Campylobacter
                fi
                #each species has its own seprate MLST
                SPEC="pubmlst_campylobacter_nonjejuni_seqdef"
                # fetus
                update_mlst_campylo "/home/external_databases/mlst/Campylobacter/fetus" "${SPEC}" "C. fetus MLST"
                # helveticus
                update_mlst_campylo "/home/external_databases/mlst/Campylobacter/helveticus" "${SPEC}" "C. helveticus MLST"
                # concisus
                update_mlst_campylo "/home/external_databases/mlst/Campylobacter/concisus" "${SPEC}" "C. concisus/curvus MLST"
                # hyointestinalis
                update_mlst_campylo "/home/external_databases/mlst/Campylobacter/hyointestinalis" "${SPEC}" "C. hyointestinalis MLST"
                # upsaliensis
                update_mlst_campylo "/home/external_databases/mlst/Campylobacter/upsaliensis" "${SPEC}" "C. upsaliensis MLST"
                # lari
                update_mlst_campylo "/home/external_databases/mlst/Campylobacter/lari" "${SPEC}" "C. lari MLST"
                # insulaenigrae
                update_mlst_campylo "/home/external_databases/mlst/Campylobacter/insulaenigrae" "${SPEC}" "C. insulaenigrae MLST"
                # lanienae
                update_mlst_campylo "/home/external_databases/mlst/Campylobacter/lanienae" "${SPEC}" "C. lanienae MLST"
                # sputorum
                update_mlst_campylo "/home/external_databases/mlst/Campylobacter/sputorum" "${SPEC}" "C. sputorum MLST"
                SPEC="pubmlst_campylobacter_seqdef"
                # jejuni
                update_mlst_campylo "/home/external_databases/mlst/Campylobacter/jejuni" "${SPEC}" "MLST"
	fi

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
	update_amrfinder >> /dev/null 2>&1
	update_kmerfinder >> /dev/null 2>&1
	update_metaphlan >> /dev/null 2>&1
	update_cge pointfinder_db >> /dev/null 2>&1
	update_cge disinfinder_db >> /dev/null 2>&1
	update_cge mlst_db >> /dev/null 2>&1
	update_cge plasmidfinder_db >> /dev/null 2>&1
	update_cge resfinder_db >> /dev/null 2>&1
	update_cge spifinder_db >> /dev/null 2>&1
	update_cge virulencefinder_db >> /dev/null 2>&1
	update_vfdb >> /dev/null 2>&1
	update_mlst ${genus}  >> /dev/null 2>&1
elif [ ${db_name} == "kraken2" ]; then
	update_kraken2 "$kraken_type" >> /dev/null 2>&1
elif [ ${db_name} == "pangolin" ]; then
	update_pangolin >> /dev/null 2>&1
elif [ ${db_name} == "freyja" ]; then
	update_freyja >> /dev/null 2>&1
elif [ ${db_name} == "nextclade" ]; then
	update_nextclade >> /dev/null 2>&1
elif [ ${db_name} == "amrfinder_plus" ]; then
	update_amrfinder >> /dev/null 2>&1
elif [ ${db_name} == "kmerfinder" ]; then
        update_kmerfinder >> /dev/null 2>&1 
elif [ ${db_name} == "metaphlan" ]; then
        update_metaphlan >> /dev/null 2>&1 
elif [ ${db_name} == "pointfinder" ]; then
        update_cge pointfinder_db >> /dev/null 2>&1
elif [ ${db_name} == "disinfinder" ]; then
        update_cge disinfinder_db >> /dev/null 2>&1
elif [ ${db_name} == "mlstfinder" ]; then
        update_cge mlst_db >> /dev/null 2>&1
elif [ ${db_name} == "plasmidfinder" ]; then
        update_cge plasmidfinder_db >> /dev/null 2>&1
elif [ ${db_name} == "resfinder" ]; then
        update_cge resfinder_db >> /dev/null 2>&1
elif [ ${db_name} == "spifinder" ]; then
        update_cge spifinder_db >> /dev/null 2>&1
elif [ ${db_name} == "virulencefinder" ]; then
        update_cge virulencefinder_db >> /dev/null 2>&1
elif [ ${db_name} == "vfdb" ]; then
	update_vfdb  >> /dev/null 2>&1 
elif [ ${db_name} == "mlst" ]; then
        update_mlst ${genus} >> /dev/null 2>&1 
fi
