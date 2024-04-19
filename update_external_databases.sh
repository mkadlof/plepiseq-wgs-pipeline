#!/usr/bin/env bash

# This script should be run periodically to update the external databases used by the pipeline.
#
# Example crontab entries:

# Update nextclade every Saturday at 3:00 AM
# Update pangolin every Saturday at 3:05 AM
# Update kraken2 every 3 months on the 1st day of the month at 3:10 AM

# 0 3 * * 6 cd /path/to/sars-illumina && bin/update_external_databases.sh nextclade
# 5 3 * * 6 cd /path/to/sars-illumina && bin/update_external_databases.sh pangolin
# 10 3 * * 6 cd /path/to/sars-illumina && bin/update_external_databases.sh freyja
# 15 3 1 */3 * cd /path/to/sars-illumina && bin/update_external_databases.sh kraken

CONTAINER="nf_illumina_sars-3.0-updater:latest"

# check if container exists
container_id=$(docker images -q $CONTAINER)
if [ -z "$container_id" ]; then
    echo "Missing container $CONTAINER. Build it before running this script!"
    exit 1
fi

# Check if argument is valid
if [ "$1" != "nextclade" ] && [ "$1" != "pangolin" ] && [ "$1" != "kraken" ] && [ "$1" != "freyja" ]
then
    echo "Invalid argument supplied. Please provide one of the following: nextclade, pangolin, kraken, freyja"
    exit 1
fi

docker run \
       --volume $(pwd)/data/${1}:/home/${1}:rw \
       --user $(id -u):$(id -g) \
       --name nf_illumina_sars-3.0-updating-${1} \
       --rm \
       $CONTAINER $1