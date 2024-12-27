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

# Function to display help message
function show_help() {
    echo "Usage: $0 --database <pangolin|nextclade|kraken2|freyja> --output-path <path> [--kraken-type <type>]"
    echo
    echo "Options:"
    echo "  --database      Name of the database to update (pangolin, nextclade, kraken2, freyja)."
    echo "  --output-path   Path to save the downloaded database."
    echo "  --kraken-type   Type of Kraken database (required if database is kraken)."
    echo "                 Possible values: standard, standard_08gb, standard_16gb, viral, minusb,"
    echo "                 pluspf, pluspf_08gb, pluspf_16gb, pluspfp, pluspfp_08gb, pluspfp_16gb,"
    echo "                 nt, eupathdb48."
    exit 1
}

# Parse arguments
DATABASE=""
OUTPUT_PATH=""
KRAKEN_TYPE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --database)
            DATABASE="$2"
            shift 2
            ;;
        --output-path)
            OUTPUT_PATH="$2"
            shift 2
            ;;
        --kraken-type)
            KRAKEN_TYPE="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            ;;
        *)
            echo "Unknown argument: $1"
            show_help
            ;;
    esac
done

# Validate required arguments
if [[ -z "$DATABASE" ]]; then
    echo "Error: --database is required."
    show_help
fi

if [[ -z "$OUTPUT_PATH" ]]; then
    echo "Error: --output-path is required."
    show_help
fi

# Validate Kraken type if database is kraken
if [[ "$DATABASE" == "kraken2" && -z "$KRAKEN_TYPE" ]]; then
    echo "Error: --kraken-type is required when --database is 'kraken2'."
    show_help
fi

if [[ "$DATABASE" == "kraken2" ]]; then
    VALID_TYPES=("standard" "standard_08gb" "standard_16gb" "viral" "minusb" \
                 "pluspf" "pluspf_08gb" "pluspf_16gb" "pluspfp" "pluspfp_08gb" \
                 "pluspfp_16gb" "nt" "eupathdb48")
    if [[ ! " ${VALID_TYPES[@]} " =~ " ${KRAKEN_TYPE} " ]]; then
        echo "Error: Invalid --kraken-type value."
        show_help
    fi
fi

# Output the parsed arguments
echo "Database: $DATABASE"
echo "Output Path: $OUTPUT_PATH"
[[ "$DATABASE" == "kraken2" ]] && echo "Kraken Type: $KRAKEN_TYPE"

# Set default path if not provided
if [ -z "$2" ]; then
    PATH_TO_USE="${OUTPUT_PATH}/${DATABASE}"
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

PATH_TO_USE=$(realpath ${PATH_TO_USE})

CONTAINER="nf_illumina_sars-4.1-updater:latest"

docker run \
       --volume ${PATH_TO_USE}:/home/external_databases:rw \
       --user $(id -u):$(id -g) \
       --name nf_illumina_sars-3.0-updating-${1} \
       --rm \
       $CONTAINER ${DATABASE} ${KRAKEN_TYPE}