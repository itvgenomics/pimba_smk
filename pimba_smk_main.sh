#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: $0 -p <PIMBA prepare mode> -r <PIMBA run mode> -g <PIMBA plot mode>"
    echo "  -p PIMBA prepare mode: paired_end, single_index, dual_index, or no"
    echo "  -r PIMBA run mode: specify the database name"
    echo "  -g PIMBA plot mode: yes or no"
    exit 1
}

# Parse command-line arguments
while getopts ":p:r:g:" opt; do
    case ${opt} in
        p )
            prepare_mode=$OPTARG
            ;;
        r )
            run_mode=$OPTARG
            ;;
        g )
            plot_mode=$OPTARG
            ;;
        \? )
            echo "Invalid option: -$OPTARG" 1>&2
            show_help
            ;;
        : )
            echo "Invalid option: -$OPTARG requires an argument" 1>&2
            show_help
            ;;
    esac
done

# Check if all arguments are provided
if [ -z "$prepare_mode" ] || [ -z "$run_mode" ] || [ -z "$plot_mode" ]; then
    echo "All arguments must be provided"
    show_help
fi

# Check PIMBA prepare mode
if [[ "$prepare_mode" != "paired_end" && "$prepare_mode" != "single_index" && "$prepare_mode" != "dual_index" && "$prepare_mode" != "no" ]]; then
    echo "Invalid PIMBA prepare mode: $prepare_mode"
    show_help
fi

# Check PIMBA run mode
if [ -z "$run_mode" ]; then
    echo "PIMBA run mode (database name) must be specified"
    show_help
fi

# Check PIMBA plot mode
if [[ "$plot_mode" != "yes" && "$plot_mode" != "no" ]]; then
    echo "Invalid PIMBA plot mode: $plot_mode"
    show_help
fi

# Implement the actions based on the arguments
if [ "$prepare_mode" == "paired_end" ]; then
    echo "Running PIMBA in paired_end prepare mode"
    snakemake --snakefile workflow/Snakefile_prepare_paired --use-singularity --configfile config/config.yaml -j 8
elif [ "$prepare_mode" == "single_index" ]; then
    echo "Running PIMBA in single_index prepare mode"
    snakemake --snakefile workflow/Snakefile_prepare_single_index --use-singularity --configfile config/config.yaml -j 8
elif [ "$prepare_mode" == "dual_index" ]; then
    echo "Running PIMBA in dual_index prepare mode"
    snakemake --snakefile workflow/Snakefile_prepare_dual_index --use-singularity --configfile config/config.yaml -j 8
elif [ "$prepare_mode" == "no" ]; then
    echo "No preparation mode selected"
fi

if [ ! -z "$run_mode" ]; then
    echo "Running PIMBA with database: $run_mode"
    snakemake --snakefile workflow/Snakefile_run --use-singularity --configfile config/config.yaml --cores 8 --singularity-args "-B /home/tfleao/Desktop/ITV/pimba_training/pimba/taxdump/:/taxdump"
fi

if [ "$plot_mode" == "yes" ]; then
    echo "Plotting graphs as specified"
    snakemake --snakefile workflow/Snakefile_plot --use-singularity --configfile config/config.yaml --cores 8
elif [ "$plot_mode" == "no" ]; then
    echo "Skipping graph plotting as specified"
fi

echo "Script execution completed"
