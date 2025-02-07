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
		rm -rf /home/external_databases/${db}/
	fi
	cd /home/external_databases
	git clone https://bitbucket.org/genomicepidemiology/${db}/
	cd ${db}
	python3 INSTALL.py /home/kma/kma_index non_interactive >> log 2>&1
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
        python3 /home/update/download_vfdb.py  >> log 2>&1
}

# MLST data
## No update mechanism
## Data for different genuses originate from eirther pubmlst or enterobase
## For Campylobacter for each species we need to call the script separetly
update_mlst_campylo() {
# subfunction for update_mlst
	local directory=${1}
	local spec=${2}
	local db=${3}
	if [ -d ${directory} ]; then
		rm -f ${directory}/*
	else
		mkdir -p ${directory}
	fi
	
	cd ${directory}
	python3 /home/update/download_mlst_campylobacter.py "${spec}" "${db}" >> log 2>&1
}

update_mlst() {
	local genus=${1}
	# Salmonella Escherichia Campylobacter
	if [ ${genus} == "Campylobacter" ]; then
		if [ ! -d "/home/external_databases/mlst/Campylobacter" ]; then
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
                        rm -f /home/external_databases/mlst/Salmonella/*
                else
                        mkdir -p /home/external_databases/mlst/Salmonella
                fi
		cd /home/external_databases/mlst/Salmonella
		python3 /home/update/download_mlst_salmonella.py >> log 2>&1
	elif [ ${genus} == "Escherichia" ]; then
		if [ -d "/home/external_databases/mlst/Escherichia" ]; then
                        rm -f /home/external_databases/mlst/Escherichia/*
                else
                        mkdir -p /home/external_databases/mlst/Escherichia
                fi
                cd /home/external_databases/mlst/Escherichia
                python3 /home/update/download_mlst_escherichia.py >> log 2>&1
	elif [ ${genus} == "all" ]; then
                if [ -d "/home/external_databases/mlst/Salmonella" ]; then
                        rm -f /home/external_databases/mlst/Salmonella/*
                else
                        mkdir -p /home/external_databases/mlst/Salmonella
                fi
                cd /home/external_databases/mlst/Salmonella
                python3 /home/update/download_mlst_salmonella.py >> log 2>&1

		if [ -d "/home/external_databases/mlst/Escherichia" ]; then
                        rm -f /home/external_databases/mlst/Escherichia/*
                else
                        mkdir -p /home/external_databases/mlst/Escherichia
                fi
                cd /home/external_databases/mlst/Escherichia
                python3 /home/update/download_mlst_escherichia.py >> log 2>&1

		if [ ! -d "/home/external_databases/mlst/Campylobacter" ]; then
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

# cgMLST related data
## like MLST data come from different sources
## for campylobacter only jejuni has cgMLST and all vrequired ariables are hardcoded within a scipt

update_cgmlst() {
	local genus=$1
	local cpus=50 # use that many workers used only for Salmonella/Escherichia data
	if [ ${genus} == "Campylobacter" ]; then
		if [ -d "/home/external_databases/cgmlst/Campylobacter/jejuni" ]; then
                        rm -f /home/external_databases/cgmlst/Campylobacter/jejuni/*
                else
                        mkdir -p /home/external_databases/cgmlst/Campylobacter/jejuni/
                fi
		cd /home/external_databases/cgmlst/Campylobacter/jejuni/
		python3 /home/update/download_cgmlst_pubmlst.py ${cpus} >> log 2>&1
	elif [ ${genus} == "Salmonella" ]; then
                if [ -d "/home/external_databases/cgmlst/Salmonella" ]; then
                        rm -f /home/external_databases/cgmlst/Salmonella/*
                else
                        mkdir -p /home/external_databases/cgmlst/Salmonella
                fi
		cd /home/external_databases/cgmlst/Salmonella
		DATABASE="senterica"
		scheme_name="cgMLST_v2"
		scheme_dir="Salmonella.cgMLSTv2"
		python3 /home/update/download_cgmlst_enterobase.py "$DATABASE" "${scheme_name}" "${scheme_dir}" ${cpus} >> log 2>&1
        elif [ ${genus} == "Escherichia" ]; then
                if [ -d "/home/external_databases/cgmlst/Escherichia" ]; then
                        rm -f /home/external_databases/cgmlst/Escherichia/*
                else
                        mkdir -p /home/external_databases/cgmlst/Escherichia
                fi
		cd /home/external_databases/cgmlst/Escherichia
		DATABASE="ecoli"
		scheme_name="cgMLST" 
		scheme_dir="Escherichia.cgMLSTv1"
		python3 /home/update/download_cgmlst_enterobase.py "$DATABASE" "${scheme_name}" "${scheme_dir}" ${cpus} >> log 2>&1
        elif [ ${genus} == "all" ]; then
		echo i"Downloading data for Escherichia at: $(date +"%H:%M %d-%m-%Y")" >> log
		if [ -d "/home/external_databases/cgmlst/Escherichia" ]; then
                        rm -f /home/external_databases/cgmlst/Escherichia/*
                else
                        mkdir -p /home/external_databases/cgmlst/Escherichia
                fi
                cd /home/external_databases/cgmlst/Escherichia
                DATABASE="ecoli"
                scheme_name="cgMLST"
                scheme_dir="Escherichia.cgMLSTv1"
                python3 /home/update/download_cgmlst_enterobase.py "$DATABASE" "${scheme_name}" "${scheme_dir}" ${cpus} >> log 2>&1
		
		echo "Downloading data for Salmonella at: $(date +"%H:%M %d-%m-%Y")" >> log
		if [ -d "/home/external_databases/cgmlst/Salmonella" ]; then
                        rm -f /home/external_databases/cgmlst/Salmonella/*
                else
                        mkdir -p /home/external_databases/cgmlst/Salmonella
                fi
                cd /home/external_databases/cgmlst/Salmonella
                DATABASE="senterica"
                scheme_name="cgMLST_v2"
                scheme_dir="Salmonella.cgMLSTv2"
                python3 /home/update/download_cgmlst_enterobase.py "$DATABASE" "${scheme_name}" "${scheme_dir}" ${cpus} >> log 2>&1

		echo "Downloading data for Campylobacter at: $(date +"%H:%M %d-%m-%Y")" >> log
		if [ -d "/home/external_databases/cgmlst/Campylobacter/jejuni" ]; then
                        rm -f /home/external_databases/cgmlst/Campylobacter/jejuni/*
                else
                        mkdir -p /home/external_databases/cgmlst/Campylobacter/jejuni/
                fi
                cd /home/external_databases/cgmlst/Campylobacter/jejuni/
                python3 /home/update/download_cgmlst_pubmlst.py ${cpus} >> log 2>&1
	fi
}

# Downloading data regarding known strains from enterobase
## There is an update mechanism so we do not remove files
update_enterobase() {
	local genus=${1}

	if [ ${genus} == "Escherichia" ]; then
                if [ ! -d "/home/external_databases/enterobase/Escherichia" ]; then
			mkdir -p /home/external_databases/enterobase/Escherichia
                fi
		cd /home/external_databases/enterobase/Escherichia
		DATABASE="ecoli" 
		CGNAME="cgMLST" 
		python3 /home/update/download_enterobase_data.py "${DATABASE}" "${CGNAME}" >> log 2>&1
        elif [ ${genus} == "Salmonella" ]; then
		if [ ! -d "/home/external_databases/enterobase/Salmonella" ]; then
                        mkdir -p /home/external_databases/enterobase/Salmonella
                fi
                cd /home/external_databases/enterobase/Salmonella
		DATABASE="senterica" 
		CGNAME="cgMLST_v2" 
		python3 /home/update/download_enterobase_data.py "${DATABASE}" "${CGNAME}" >> log 2>&1
	elif [[ ${genus} == "all" || ${genus} == "Campylobacter" ]]; then
		if [ ! -d "/home/external_databases/enterobase/Escherichia" ]; then
                        mkdir -p /home/external_databases/enterobase/Escherichia
                fi
                cd /home/external_databases/enterobase/Escherichia
                DATABASE="ecoli"
                CGNAME="cgMLST"
                python3 /home/update/download_enterobase_data.py "${DATABASE}" "${CGNAME}" >> log 2>&1

		if [ ! -d "/home/external_databases/enterobase/Salmonella" ]; then
                        mkdir -p /home/external_databases/enterobase/Salmonella
                fi
                cd /home/external_databases/enterobase/Salmonella
                DATABASE="senterica"
                CGNAME="cgMLST_v2"
                python3 /home/update/download_enterobase_data.py "${DATABASE}" "${CGNAME}" >> log 2>&1
	fi

}

# Downloading data to get analogue of pHierCC for C.jejuni and for historical analysis
## There is an update mechanism

update_pubmlst() {
	if [ ! -d "/home/external_databases/pubmlst/Campylobacter/jejuni" ]; then
		mkdir -p /home/external_databases/pubmlst/Campylobacter/jejuni
	fi
	cd /home/external_databases/pubmlst/Campylobacter/jejuni
	python3 /home/update/download_pubmlst_data.py >> log 2>&1

}


# Data from custom clustering of cgMLST data
## No update mechanism so we donload again files
update_phiercc() {
	local genus=${1}
	if [ ${genus} == "Campylobacter" ]; then
		if [ -d "/home/external_databases/phiercc_local/Campylobacter/jejuni" ]; then
			rm -rf /home/external_databases/phiercc_local/Campylobacter/jejuni/*
		else
			mkdir -p /home/external_databases/phiercc_local/Campylobacter/jejuni/
		fi
		cd /home/external_databases/phiercc_local/Campylobacter/jejuni/
                
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Campylobacter/profile_complete_linkage.HierCC.gz
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Campylobacter/profile_complete_linkage.HierCC.index
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Campylobacter/profile_single_linkage.HierCC.gz
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Campylobacter/profile_single_linkage.HierCC.index

	elif [ ${genus} == "Escherichia" ]; then
		if [ -d "/home/external_databases/phiercc_local/Escherichia" ]; then
                        rm -rf /home/external_databases/phiercc_local/Escherichia/*
		else
			mkdir -p /home/external_databases/phiercc_local/Escherichia
		fi
		cd /home/external_databases/phiercc_local/Escherichia
		
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Escherichia/profile_single_linkage.HierCC.gz
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Escherichia/profile_single_linkage.HierCC.index
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Escherichia/profile_complete_linkage.HierCC.gz
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Escherichia/profile_complete_linkage.HierCC.index

	elif [ ${genus} == "Salmonella" ]; then
		if [ -d "/home/external_databases/phiercc_local/Salmonella" ]; then
			rm -rf /home/external_databases/phiercc_local/Salmonella/*
		else
			mkdir -p /home/external_databases/phiercc_local/Salmonella
		fi
		cd /home/external_databases/phiercc_local/Salmonella
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Salmonella/profile_complete_linkage.HierCC.gz
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Salmonella/profile_complete_linkage.HierCC.index
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Salmonella/profile_single_linkage.HierCC.gz
		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Salmonella/profile_single_linkage.HierCC.index
	
	
	elif [ ${genus} == "all" ]; then
		if [ -d "/home/external_databases/phiercc_local/Salmonella" ]; then
                        rm -rf /home/external_databases/phiercc_local/Salmonella/*
                else
                        mkdir -p /home/external_databases/phiercc_local/Salmonella
                fi

		cd /home/external_databases/phiercc_local/Salmonella

		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Salmonella/profile_complete_linkage.HierCC.gz
                wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Salmonella/profile_complete_linkage.HierCC.index
                wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Salmonella/profile_single_linkage.HierCC.gz
                wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Salmonella/profile_single_linkage.HierCC.index

		if [ -d "/home/external_databases/phiercc_local/Escherichia" ]; then
                        rm -rf /home/external_databases/phiercc_local/Escherichia/*
                else
                        mkdir -p /home/external_databases/phiercc_local/Escherichia
                fi
                cd /home/external_databases/phiercc_local/Escherichia

		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Escherichia/profile_single_linkage.HierCC.gz
                wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Escherichia/profile_single_linkage.HierCC.index
                wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Escherichia/profile_complete_linkage.HierCC.gz
                wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Escherichia/profile_complete_linkage.HierCC.index

		if [ -d "/home/external_databases/phiercc_local/Campylobacter/jejuni" ]; then
			rm -rf /home/external_databases/phiercc_local/Campylobacter/jejuni/*
		else
			mkdir -p /home/external_databases/phiercc_local/Campylobacter/jejuni/
		fi
		cd /home/external_databases/phiercc_local/Campylobacter/jejuni/

		wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Campylobacter/profile_complete_linkage.HierCC.gz
                wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Campylobacter/profile_complete_linkage.HierCC.index
                wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Campylobacter/profile_single_linkage.HierCC.gz
                wget https://github.com/michallaz/phiercc_pzh_mod/raw/refs/heads/main/plepiseq_data/Campylobacter/profile_single_linkage.HierCC.index

	fi
}

# Uniref50 and Uniref50 for Virual sequences used by alphafold
## No update mechanism
update_alphafold() {
	if [ -d "/home/external_databases/alphafold/uniref_viruses" ]; then
		rm -rf /home/external_databases/alphafold/uniref50/*
	else
		mkdir -p /home/external_databases/alphafold/uniref50/
	fi
	
	wget -O /home/external_databases/alphafold/uniref50/uniref50_viral.fasta "https://rest.uniprot.org/uniref/stream?format=fasta&query=%28%28taxonomy_id%3A10239%29+AND+%28identity%3A0.5%29%29"
	wget -O /home/external_databases/alphafold/uniref50/uniref50.fasta.gz https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
	gunzip /home/external_databases/alphafold/uniref50/uniref50.fasta.gz

	if [ -d "/home/external_databases/alphafold/uniprot" ]; then
                rm -rf /home/external_databases/alphafold/uniprot/uniprot_sprot.fasta
        else
                mkdir -p /home/external_databases/alphafold/uniprot
        fi

	wget -O /home/external_databases/alphafold/uniprot/uniprot_sprot.fasta.gz "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
	gunzip /home/external_databases/alphafold/uniprot/uniprot_sprot.fasta.gz
}
#############
# Main code *
#############


db_name=$1
kraken_type=$2
genus=$3
if [ ${db_name} == "all" ];then
        echo "Downloading data for kraken2 at: $(date +"%H:%M %d-%m-%Y")"
	update_kraken2 "$kraken_type" >> /dev/null 2>&1
	echo "Downloading data for pangolin at: $(date +"%H:%M %d-%m-%Y")"
        update_pangolin >> /dev/null 2>&1
	echo "Downloading data for freyja at: $(date +"%H:%M %d-%m-%Y")"
        update_freyja >> /dev/null 2>&1
	echo "Downloading data for nextclade at: $(date +"%H:%M %d-%m-%Y")"
        update_nextclade >> /dev/null 2>&1
	echo "Downloading data for AMRfinder_plus at: $(date +"%H:%M %d-%m-%Y")"
	update_amrfinder >> /dev/null 2>&1
	echo "Downloading data for kmerfinder at: $(date +"%H:%M %d-%m-%Y")"
	update_kmerfinder >> /dev/null 2>&1
	echo "Downloading data for metaphlan at: $(date +"%H:%M %d-%m-%Y")"
	update_metaphlan >> /dev/null 2>&1
	echo "Downloading data for pointfinder at: $(date +"%H:%M %d-%m-%Y")"
	update_cge pointfinder_db >> /dev/null 2>&1
	echo "Downloading data for disinfinder at: $(date +"%H:%M %d-%m-%Y")"
	update_cge disinfinder_db >> /dev/null 2>&1
	echo "Downloading data for mlst_db at: $(date +"%H:%M %d-%m-%Y")"
	update_cge mlst_db >> /dev/null 2>&1
	echo "Downloading data for plasmidfinder at: $(date +"%H:%M %d-%m-%Y")"
	update_cge plasmidfinder_db >> /dev/null 2>&1
	echo "Downloading data for resfinder at: $(date +"%H:%M %d-%m-%Y")"
	update_cge resfinder_db >> /dev/null 2>&1
	echo "Downloading data for spifinder at: $(date +"%H:%M %d-%m-%Y")"
	update_cge spifinder_db >> /dev/null 2>&1
	echo "Downloading data for virulencefinder at: $(date +"%H:%M %d-%m-%Y")"
	update_cge virulencefinder_db >> /dev/null 2>&1
	echo "Downloading data for vfcb at: $(date +"%H:%M %d-%m-%Y")"
	update_vfdb >> /dev/null 2>&1
	echo "Downloading MLST data at: $(date +"%H:%M %d-%m-%Y")"
	update_mlst ${genus}  >> /dev/null 2>&1
	echo "Downloading cgMLST data at: $(date +"%H:%M %d-%m-%Y")"
	update_cgmlst ${genus}  >> /dev/null 2>&1
	echo "Downloading pubmlst data at: $(date +"%H:%M %d-%m-%Y")"
	update_pubmlst >> /dev/null 2>&1
	echo "Downloading enterobase data at: $(date +"%H:%M %d-%m-%Y")"
	update_enterobase ${genus} >> /dev/null 2>&1
	echo "Downloading hiercc data at: $(date +"%H:%M %d-%m-%Y")"
	update_phiercc ${genus} >> /dev/null 2>&1
	echo "Downloading swissprot and uniref50 for viral sequences at: $(date +"%H:%M %d-%m-%Y")"
	update_alphafold >> /dev/null 2>&1
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
elif [ ${db_name} == "cgmlst" ]; then
	update_cgmlst ${genus} >> /dev/null 2>&1
elif [ ${db_name} == "pubmlst" ]; then
	update_pubmlst >> /dev/null 2>&1
elif [ ${db_name} == "enterobase" ]; then
	update_enterobase ${genus}
elif [ ${db_name} == "phiercc" ]; then
	update_phiercc ${genus} >> /dev/null 2>&1
elif [ ${db_name} == "alphafold" ]; then
	update_alphafold >> /dev/null 2>&1
fi
