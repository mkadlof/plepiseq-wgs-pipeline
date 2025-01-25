#!/bin/bash

# Skrypt sciaga aktualna baze z serwera http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/
# Uwaga baza MetaPhlana jest duza > 60Gb 

# This path is HARDCODED
cd /home/external_databases/metaphlan

if [ -e mpa_latest ]; then
	OLD_PREFIX=`cat mpa_latest`
	rm mpa_latest
fi
curl -u anonymous:anonymous http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_latest -O


PREFIX=`cat mpa_latest`
if [ ${PREFIX} == ${OLD_PREFIX} ]; then
	echo "Database is up to date"
	exit 1
else
	rm mpa_v*
	rm files*
	curl -u anonymous:anonymous http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/ | grep ${PREFIX} | cut -d " " -f8 |  cut -d '"' -f2 >> files_main.txt
	curl -u anonymous:anonymous http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/ | grep ${PREFIX} | cut -d " " -f8 | cut -d '"' -f2 >> files_bowtie_indexes.txt

	for K in `cat files_main.txt`
	do
	curl -u  anonymous:anonymous http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/${K} -O
	if [ 'tar' == `echo ${K} | cut -d '.' -f2` ]; then
		tar -xf ${K}
	fi

	done

	for K in `cat files_bowtie_indexes.txt`
	do
	curl -u  anonymous:anonymous http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/${K} -O
	if [ 'tar' == `echo ${K} | cut -d '.' -f2` ]; then
        	tar -xf ${K}
	fi
	done

	# unzip the files
	for K in `ls *bz2`; do bunzip2 ${K}; done
fi
