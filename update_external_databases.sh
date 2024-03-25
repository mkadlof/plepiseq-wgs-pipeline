#!/usr/bin/env bash

# This script should be run periodically to update the external databases used by the pipeline.
#

CONTAINER="nf_illumina_sars-3.0-updater:latest"

# check if container exists
container_id=$(docker images -q $CONTAINER)
if [ -z "$container_id" ]; then
    echo "Missing container $CONTAINER. Build it before running this script!"
    exit 1
fi

docker run \
       --volume $(pwd)/data/nextclade:/home/nextclade \
       --volume $(pwd)/data/pangolin:/home/pangolin \
       --user $(id -u):$(id -g) \
       $CONTAINER