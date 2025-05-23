#!/bin/bash

SECONDS=0

# Function to display help message
show_help() {
    echo "Usage: $0 -p <PIMBA prepare mode> -r <PIMBA run mode> -g <PIMBA plot mode> -t <number of threads> -c <config file> -d <work directory>"
    echo "  -p PIMBA prepare mode: paired_end, single_index, dual_index, or no"
    echo "  -r PIMBA run mode: specify the database name"
    echo "  -g PIMBA plot mode: yes or no"
    echo "  -t number of threads: specify the number of threads to use"
    echo "  -c config file: specify the path to the config file"
    echo "  -d work directory: specify the working directory"
    exit 1
}

# Parse command-line arguments
while getopts ":p:r:g:t:c:d:l:u" opt; do
    case ${opt} in
        p ) prepare_mode=$OPTARG ;;
        r ) run_mode=$OPTARG ;;
        g ) plot_mode=$OPTARG ;;
        t ) threads=$OPTARG ;;
        c ) config_file=$OPTARG ;;
        d ) workdir=$OPTARG ;;
        l ) place_mode=$OPTARG ;;
        u ) SETUNLOCK="--unlock" ;;
        \? ) echo "Invalid option: -$OPTARG" 1>&2; show_help ;;
        : ) echo "Invalid option: -$OPTARG requires an argument" 1>&2; show_help ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$prepare_mode" ] || [ -z "$run_mode" ] || [ -z "$plot_mode" ] || [ -z "$threads" ] || [ -z "$config_file" ] || [ -z "$workdir" ] || [ -z "$place_mode" ]; then
    echo "All arguments must be provided"
    show_help
fi

# Convert file path to absolute path
config_file=$(realpath "$config_file")
workdir=$(realpath "$workdir")

# Extract necessary paths from the config file
rawdatadir=$(grep '^rawdatadir:' "$config_file" | awk '{print $2}' | tr -d "'")
raw_fastq_single=$(grep '^raw_fastq_single:' "$config_file" | awk '{print $2}' | tr -d "'")
raw_fastq_dual=$(grep '^raw_fastq_dual:' "$config_file" | awk '{print $2}' | tr -d "'")
ncbi_db=$(grep '^NCBI-DB:' "$config_file" | awk '{print $2}' | tr -d "'")
taxdump=$(grep '^taxdump:' "$config_file" | awk '{print $2}' | tr -d "'")
remote=$(grep '^remote:' "$config_file" | awk '{print $2}' | tr -d "'")
metadata=$(grep '^metadata:' "$config_file" | awk '{print $2}' | tr -d "'")
adapters=$(grep '^adapters:' "$config_file" | awk '{print $2}' | tr -d "'")

# ---------------- PIMBA Prepare ----------------
if [ "$prepare_mode" != "no" ]; then
    if [ "$prepare_mode" == "paired_end" ]; then
        if [ -z "$rawdatadir" ] || [ ! -d "$rawdatadir" ]; then
            echo "ERROR: rawdatadir is not defined or does not exist."
            exit 1
        fi
        if [ -z "$adapters" ] || [ ! -e "$adapters" ]; then
            echo "ERROR: adapters file is not defined or does not exist."
            exit 1
        fi
        echo "Running PIMBA in paired_end prepare mode"
        snakemake --snakefile workflow/Snakefile_prepare_paired --use-singularity \
            --configfile "$config_file" --cores "$threads" --directory "$workdir" \
            --singularity-args "-B $rawdatadir -B $adapters" ${SETUNLOCK}
    
    elif [ "$prepare_mode" == "single_index" ]; then
        if [ -z "$raw_fastq_single" ] || [ ! -e "$raw_fastq_single" ]; then
            echo "ERROR: raw_fastq_single is not defined or does not exist."
            exit 1
        fi
        raw_fastq_dir=$(dirname "$raw_fastq_single")
        echo "Running PIMBA in single_index prepare mode"
        snakemake --snakefile workflow/Snakefile_prepare_single_index --use-singularity \
            --configfile "$config_file" --cores "$threads" --directory "$workdir" \
            --singularity-args "-B $raw_fastq_dir" ${SETUNLOCK}

    elif [ "$prepare_mode" == "dual_index" ]; then
        if [ -z "$raw_fastq_dual" ] || [ ! -e "$raw_fastq_dual" ]; then
            echo "ERROR: raw_fastq_dual is not defined or does not exist."
            exit 1
        fi
        raw_fastq_dir=$(dirname "$raw_fastq_dual")
        echo "Running PIMBA in dual_index prepare mode"
        snakemake --snakefile workflow/Snakefile_prepare_dual_index --use-singularity \
            --configfile "$config_file" --cores "$threads" --directory "$workdir" \
            --singularity-args "-B $raw_fastq_dir" ${SETUNLOCK}
    fi
else
    echo "Skipping PIMBA Prepare as specified"
fi

# ---------------- PIMBA Run ----------------
if [ "$run_mode" != "no" ]; then
    case "$run_mode" in
        "16S-NCBI"|"COI-NCBI"|"ITS-PLANTS-NCBI"|"ITS-FUNGI-NCBI"|"ALL-NCBI")
            if [ "$remote" == "yes" ]; then
                if [ -z "$taxdump" ] || [ ! -d "$taxdump" ]; then
                    echo "ERROR: taxdump is not defined or does not exist."
                    exit 1
                fi
                echo "Running PIMBA with database $run_mode plus remote mode as $remote"
                snakemake --snakefile workflow/Snakefile_run --use-singularity \
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
                echo "Running PIMBA with database $run_mode plus remote mode as $remote"
                snakemake --snakefile workflow/Snakefile_run --use-singularity \
                    --configfile "$config_file" --cores "$threads" --directory "$workdir" \
                    --singularity-args "-B $ncbi_db -B $taxdump:/taxdump -B $taxdump" ${SETUNLOCK}
            fi
            ;;
        *)
            db_path=$(grep "^$run_mode-DB:" "$config_file" | awk '{print $2}' | tr -d "'")
            if [ -z "$db_path" ]; then
                echo "Database path for $run_mode-DB not found in config file. Using custom database"
                db_path=$run_mode
            fi
            if [ ! -d "$db_path" ]; then
                echo "ERROR: Database path '$db_path' does not exist."
                exit 1
            fi
            echo "Running PIMBA with database: $db_path"
            snakemake --snakefile workflow/Snakefile_run --use-singularity \
                --configfile "$config_file" --cores "$threads" --directory "$workdir" \
                --singularity-args "-B $db_path" ${SETUNLOCK}
            ;;
    esac
else
    echo "Skipping PIMBA Run as specified"
fi

# ---------------- PIMBA Plot ----------------
if [ "$plot_mode" == "yes" ]; then
    if [ -z "$metadata" ] || [ ! -e "$metadata" ]; then
        echo "ERROR: metadata is not defined or does not exist."
        exit 1
    fi
    echo "Plotting graphs as specified"
    snakemake --snakefile workflow/Snakefile_plot --use-singularity \
        --configfile "$config_file" --cores "$threads" --directory "$workdir" \
        --singularity-args "-B $metadata" ${SETUNLOCK}
else
    echo "Skipping graph plotting as specified"
fi

# ---------------- PIMBA Place ----------------

if [ "$place_mode" == "yes" ]; then
    snakemake --snakefile workflow/Snakefile_place --cores "$threads" --use-conda --configfile config/config_place.yaml --conda-frontend conda ${SETUNLOCK}
else
    echo "Skipping PIMBA Place as specified"
fi

duration=$SECONDS
echo "Script execution completed in $((duration / 60)) minutes and $((duration % 60)) seconds"