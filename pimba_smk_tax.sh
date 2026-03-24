#!/bin/bash

SECONDS=0

# Function to display help message
show_help() {
    echo "Usage: $0 -r <PIMBA tax mode> -t <threads> -c <config file> -d <work directory>"
    echo "  -r PIMBA tax mode: specify the database name"
    echo "  -t number of threads"
    echo "  -c config file"
    echo "  -d work directory"
    exit 1
}

# Parse command-line arguments
while getopts ":r:t:c:d:u" opt; do
    case ${opt} in
        r ) tax_mode=$OPTARG ;;
        t ) threads=$OPTARG ;;
        c ) config_file=$OPTARG ;;
        d ) workdir=$OPTARG ;;
        u ) SETUNLOCK="--unlock" ;;
        \? ) echo "Invalid option: -$OPTARG" 1>&2; show_help ;;
        : ) echo "Option -$OPTARG requires an argument" 1>&2; show_help ;;
    esac
done

# Check required arguments
if [ -z "$tax_mode" ] || [ -z "$threads" ] || [ -z "$config_file" ] || [ -z "$workdir" ]; then
    echo "Missing required arguments"
    show_help
fi

# Convert to absolute paths
config_file=$(realpath "$config_file")
workdir=$(realpath "$workdir")

# Fetch Singularity images
python "$workdir/workflow/scripts/singularity.py" --config "$config_file"

# Extract config values
ncbi_db=$(grep '^NCBI-DB:' "$config_file" | awk '{print $2}' | tr -d "'")
taxdump=$(grep '^taxdump:' "$config_file" | awk '{print $2}' | tr -d "'")
remote=$(grep '^remote:' "$config_file" | awk '{print $2}' | tr -d "'")

# ---------------- PIMBA Tax ----------------
case "$tax_mode" in
    "16S-NCBI"|"COI-NCBI"|"ITS-PLANTS-NCBI"|"ITS-FUNGI-NCBI"|"ALL-NCBI")
        
        if [ "$remote" == "yes" ]; then
            if [ -z "$taxdump" ] || [ ! -d "$taxdump" ]; then
                echo "ERROR: taxdump is not defined or does not exist."
                exit 1
            fi

            echo "Running PIMBA with $tax_mode (remote mode = yes)"
            snakemake --snakefile workflow/Snakefile_tax --use-singularity \
                --configfile "$config_file" --cores "$threads" --directory "$workdir" \
                --singularity-args "-B $taxdump:/taxdump -B $taxdump" ${SETUNLOCK}

        else
            if [ -z "$ncbi_db" ] || [ ! -d "$ncbi_db" ]; then
                echo "ERROR: ncbi_db is not defined or does not exist."
                exit 1
            fi

            if [ -z "$taxdump" ] || [ ! -d "$taxdump" ]; then
                echo "ERROR: taxdump is not defined or does not exist."
                exit 1
            fi

            echo "Running PIMBA with $tax_mode (remote mode = no)"
            snakemake --snakefile workflow/Snakefile_tax --use-singularity \
                --configfile "$config_file" --cores "$threads" --directory "$workdir" \
                --singularity-args "-B $ncbi_db -B $taxdump:/taxdump -B $taxdump" ${SETUNLOCK}
        fi
        ;;

    *)
        db_path=$(grep "^$tax_mode-DB:" "$config_file" | awk '{print $2}' | tr -d "'")

        if [ -z "$db_path" ]; then
            echo "Database path for $run_mode-DB not found in config file. Using custom database"
            db_path=$tax_mode
        fi

        if [ ! -d "$db_path" ]; then
            echo "ERROR: Database path '$db_path' does not exist."
            exit 1
        fi

        echo "Running PIMBA with database: $db_path"
        snakemake --snakefile workflow/Snakefile_tax --use-singularity \
            --configfile "$config_file" --cores "$threads" --directory "$workdir" \
            --singularity-args "-B $db_path" ${SETUNLOCK}
        ;;
esac

duration=$SECONDS
echo "Completed in $((duration / 60)) min $((duration % 60)) sec"