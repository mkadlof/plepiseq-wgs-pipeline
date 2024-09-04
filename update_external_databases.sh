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

usage() {
    echo "Usage: $0 all|pangolin|nextclade|kraken|freyja [PATH]"
    echo "  all       - Update all databases. Default path: ./external_databases"
    echo "  pangolin  - Update the pangolin database. Default path: ./external_databases/pangolin"
    echo "  nextclade - Update the nextclade database. Default path: ./external_databases/nextclade"
    echo "  kraken    - Update the kraken database. Default path: ./external_databases/kraken"
    echo "  freyja    - Update the freyja database. Default path: ./external_databases/freyja"
    echo "  PATH      - Optional path to the database directory"
}

CONTAINER="nf_illumina_sars-3.0-updater:latest"

# Check if container exists
container_id=$(docker images -q $CONTAINER)
if [ -z "$container_id" ]; then
    echo "Missing container $CONTAINER. Build it before running this script!"
    exit 1
fi

# Check if argument is valid
if [ "$1" != "nextclade" ] && [ "$1" != "pangolin" ] && [ "$1" != "kraken" ] && [ "$1" != "freyja" ] && [ "$1" != "all" ]; then
    echo "Invalid argument supplied. Please provide one of the following: all, nextclade, pangolin, kraken, freyja"
    usage
    exit 1
fi

# Set default path if not provided
if [ -z "$2" ]; then
    if [ "$1" == "all" ]; then
        PATH_TO_USE="./external_databases"
    else
        PATH_TO_USE="./external_databases/$1"
    fi
else
    PATH_TO_USE="$2"
fi

# Ensure the directory exists
if [ ! -d "$PATH_TO_USE" ]; then
    echo "Directory $PATH_TO_USE does not exist. Creating it..."
    mkdir -p "$PATH_TO_USE"
fi

# Check if the current user has write permissions to the directory
if [ ! -w "$PATH_TO_USE" ]; then
    echo "Current user does not have write permissions to the directory $PATH_TO_USE."
    exit 1
fi

docker run \
       --volume $(pwd)/${PATH_TO_USE}:/home/external_databases:rw \
       --user $(id -u):$(id -g) \
       --name nf_illumina_sars-3.0-updating-${1} \
       --rm \
       $CONTAINER $1