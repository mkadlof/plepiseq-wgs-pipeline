#!/bin/bash

echo "CGE ftp is often busy, hence it is sometimes required to run this script several times before the data can be downloaded sucesfully"

if [ -d kmerfinder_db ]; then
	rm -rf kmerfinder_db
fi

git clone https://bitbucket.org/genomicepidemiology/kmerfinder_db.git
cd kmerfinder_db/
bash INSTALL.sh `pwd` bacteria
