#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: $0 -p <PIMBA prepare mode> -r <PIMBA run mode> -g <PIMBA plot mode> -t <number of threads> -c <config file>"
    echo "  -p PIMBA prepare mode: paired_end, single_index, dual_index, or no"
    echo "  -r PIMBA run mode: specify the database name"
    echo "  -g PIMBA plot mode: yes or no"
    echo "  -t number of threads: specify the number of threads to use"
    echo "  -c config file: specify the path to the config file"
    exit 1
}

# Parse command-line arguments
while getopts ":p:r:g:t:c:" opt; do
    case ${opt} in
        p ) prepare_mode=$OPTARG ;;
        r ) run_mode=$OPTARG ;;
        g ) plot_mode=$OPTARG ;;
        t ) threads=$OPTARG ;;
        c ) config_file=$OPTARG ;;
        \? ) echo "Invalid option: -$OPTARG" 1>&2; show_help ;;
        : ) echo "Invalid option: -$OPTARG requires an argument" 1>&2; show_help ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$prepare_mode" ] || [ -z "$run_mode" ] || [ -z "$plot_mode" ] || [ -z "$threads" ] || [ -z "$config_file" ]; then
    echo "All arguments must be provided"
    show_help
fi

# Convert config file path to absolute path
config_file=$(realpath "$config_file")

# Extract necessary paths from the config file
rawdatadir=$(grep '^rawdatadir:' "$config_file" | awk '{print $2}' | tr -d "'")
raw_fastq_single=$(grep '^raw_fastq_single:' "$config_file" | awk '{print $2}' | tr -d "'")
raw_fastq_dual=$(grep '^raw_fastq_dual:' "$config_file" | awk '{print $2}' | tr -d "'")
ncbi_db=$(grep '^NCBI-DB:' "$config_file" | awk '{print $2}' | tr -d "'")
taxdump=$(grep '^taxdump:' "$config_file" | awk '{print $2}' | tr -d "'")
remote=$(grep '^remote:' "$config_file" | awk '{print $2}' | tr -d "'")
metadata=$(grep '^metadata:' "$config_file" | awk '{print $2}' | tr -d "'")

if [ "$prepare_mode" != "no" ]; then
    if [ "$prepare_mode" == "paired_end" ]; then
        echo "Running PIMBA in paired_end prepare mode"
        snakemake --snakefile workflow/Snakefile_prepare_paired --use-singularity --configfile "$config_file" --cores "$threads" --singularity-args "-B $rawdatadir"
    elif [ "$prepare_mode" == "single_index" ]; then
        echo "Running PIMBA in single_index prepare mode"
        raw_fastq_dir=$(dirname "$raw_fastq_single")
        snakemake --snakefile workflow/Snakefile_prepare_single_index --use-singularity --configfile "$config_file" --cores "$threads" --singularity-args "-B $raw_fastq_dir"
    elif [ "$prepare_mode" == "dual_index" ]; then
        echo "Running PIMBA in dual_index prepare mode"
        raw_fastq_dir=$(dirname "$raw_fastq_dual")
        snakemake --snakefile workflow/Snakefile_prepare_dual_index --use-singularity --configfile "$config_file" --cores "$threads" --singularity-args "-B $raw_fastq_dir"
    fi
else
    echo "Skipping PIMBA prepare as specified"
fi

if [ "$run_mode" != "no" ]; then
    if [[ "$run_mode" == *"NCBI"* ]]; then
        echo "Running PIMBA with database $run_mode plus remote mode as $remote"
        if [ "$remote" == "yes" ]; then
            snakemake --snakefile workflow/Snakefile_run --use-singularity --configfile "$config_file" --cores "$threads" --singularity-args "-B $taxdump:/taxdump -B $taxdump"
        else
            snakemake --snakefile workflow/Snakefile_run --use-singularity --configfile "$config_file" --cores "$threads" --singularity-args "-B $ncbi_db -B $taxdump:/taxdump -B $taxdump"
        fi
    else
        db_path=$(grep "^$run_mode-DB:" "$config_file" | awk '{print $2}')
        if [ -z "$db_path" ]; then
            echo "Database path for $run_mode-DB not found in config file. Using custom database"
            db_path=$run_mode
        fi
        echo "Running PIMBA with database: $run_mode"
        snakemake --snakefile workflow/Snakefile_run --use-singularity --configfile "$config_file" --cores "$threads" --singularity-args "-B $db_path"
    fi
else
    echo "Skipping PIMBA run as specified"
fi

if [ "$plot_mode" == "yes" ]; then
    echo "Plotting graphs as specified"
    snakemake --snakefile workflow/Snakefile_plot --use-singularity --configfile "$config_file" --cores "$threads" --singularity-args "-B $metadata"
else
    echo "Skipping graph plotting as specified"
fi

echo "Script execution completed"