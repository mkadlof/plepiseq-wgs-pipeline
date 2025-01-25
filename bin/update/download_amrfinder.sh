#!/bin/bash
# 'SKRYPT WYMAGA BLASTA spojnego z DOCKERFILEM (aktualnie w wersji 2.13.0+)'
# download the files from an ftp server
# and index the database

# updating mechanism 
# on the ftp server there is a version.txt file with date

# This path is HARDCODED
cd /home/external_databases/amrfinder_plus

if [ -e version.txt ]; then
	# Download the data only if new version is availalbe
	OLD_VERSION=`cat version.txt`
	NEW_VERSION=`curl -u anonymous:anonymous https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/version.txt`
else
	OLD_VERSION="1990-01-01"
fi

if [ "${OLD_VERSION}" !=  "${NEW_VERSION}" ]; then
	echo "New version: ${NEW_VERSION} detected"

	curl -u anonymous:anonymous ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ | while read -r line; do
    	FILE_NAME=$(echo $line | awk '{print $9}')
    	# For now just remove the old file
    	rm ${FILE_NAME}*
    	# Download each file
    	curl -u anonymous:anonymous ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/$FILE_NAME -O
    
    	# check if file is in fasta format and whether it contains nucl or ptotein sequence
    	HEADER=`head -1 ${FILE_NAME} | grep ">" | wc -l`
    	TYPE=`head -2 ${FILE_NAME} | tail -1  | fold -w 1 | grep -v [ATGCatgc] | wc -l`
    	if [ ${TYPE} -eq 0 ]; then
		DBTYPE='nucl'
	else
		DBTYPE='prot'
	fi

    	if [ ${HEADER} -eq 1 ]; then
	    	makeblastdb -in ${FILE_NAME} -dbtype ${DBTYPE}
	fi
	done
else
	echo "No need for an update your version: ${OLD_VERSION} is up-to-date"
fi
