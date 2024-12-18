# PIMBA v3.0 in Snakemake - User Guide

Authors: Tiago Ferreira Leão, Renato R. M. Oliveira

Document Version: 1.0

Date: 13/12/2024

## Description
PIMBA (Pipeline for MetaBarcoding Analysis) (Oliveira et al. 2021) is a tool for metabarcoding analysis that allows users to create their own database, overcoming limitations of other similar tools. PIMBA adapts the Qiime/BMP pipeline (Bolyen et al. 2019)(Pylro et al. 2014) for OTU clustering and includes optional OTU corrections using the LULU algorithm (Frøslev et al. 2017). It supports paired and unpaired reads (with single or double indexing options) and enables ASV inference using Swarm (Mahé et al. 2021). PIMBA provides preliminary abundance and diversity analyses automatically. This pipeline is written in Snakemake (Köster and Rahmann 2012), a workflow management system designed to streamline the execution of complex data analysis pipelines, offering significant advantages over traditional bash scripts. Snakemake provides a structured and modular approach that enhances the readability and manageability of the pipeline. This guide will help you install and run PIMBA using Snakemake.

## References

1. Bolyen, Evan, Jai Ram Rideout, Matthew R. Dillon, Nicholas A. Bokulich, Christian C. Abnet, Gabriel A. Al-Ghalith, Harriet Alexander, et al. 2019. “Reproducible, Interactive, Scalable and Extensible Microbiome Data Science Using QIIME 2.” Nature Biotechnology 37 (8): 852–57.
2. Frøslev, Tobias Guldberg, Rasmus Kjøller, Hans Henrik Bruun, Rasmus Ejrnæs, Ane Kirstine Brunbjerg, Carlotta Pietroni, and Anders Johannes Hansen. 2017. “Algorithm for Post-Clustering Curation of DNA Amplicon Data Yields Reliable Biodiversity Estimates.” Nature Communications 8 (1): 1188.
3. Köster, Johannes, and Sven Rahmann. 2012. “Snakemake--a Scalable Bioinformatics Workflow Engine.” Bioinformatics (Oxford, England) 28 (19): 2520–22.
4. Mahé, Frédéric, Lucas Czech, Alexandros Stamatakis, Christopher Quince, Colomban de Vargas, Micah Dunthorn, and Torbjørn Rognes. 2021. “Swarm v3: Towards Tera-Scale Amplicon Clustering.” Bioinformatics (Oxford, England) 38 (1): 267–69.
5. Oliveira, Renato R. M., Raíssa Silva, Gisele L. Nunes, and Guilherme Oliveira. 2021. “PIMBA: A PIpeline for MetaBarcoding Analysis.” Advances in Bioinformatics and Computational Biology, 106–16.
6. Pylro, Victor S., Luiz Fernando W. Roesch, Daniel K. Morais, Ian M. Clark, Penny R. Hirsch, and Marcos R. Tótola. 2014. “Data Analysis for 16S Microbial Profiling from Different Benchtop Sequencing Platforms.” Journal of Microbiological Methods 107 (December):30–37.

## Prerequisites

Before installing the software, make sure you have the following:

- A Linux-based operating system (e.g., Ubuntu, CentOS, Fedora) or Windows WSL (Windows Subsystem for Linux)
- Python (version 3.5 or later) installed on your system
- Git

## Installation

### A) Anaconda

To run Snakemake, you need to install Anaconda by following these steps:

1. Download the installation file from here: https://www.anaconda.com/download/;
2. In the terminal, navigate to the directory where the installation file is located, and run:

`bash <filename>-latest-Linux-x86_64.sh`

For example:

`bash Anaconda3-2024.10.1-latest-Linux-x86_64.sh`

3. Follow the instructions on the installer screens. If you are unsure about any settings, accept the defaults. You can change them later;
4. To apply the changes, close and reopen the terminal window;
5. Test your installation. In the terminal window, run the command conda list. A list of installed packages will appear if the installation was successful.

### B) Snakemake

You also need to create a new conda environment and install Snakemake and Singularity by following these steps:

`conda create -n snakemake_env -c bioconda -c conda-forge singularity=3.8.6 snakemake=7.32.4`

`conda activate snakemake_env`

### C) Clone the GitHub Repository

Finally, you need to clone the GitHub repository containing the code to run PIMBA in Snakemake. Use the following command:

`git clone https://github.com/itvgenomics/pimba_smk.git`

`cd pimba_smk`

You should observe the following files within the cloned directory:
- Config: Directory containing the file config.yaml where all analysis parameters will be configured. For a standard user, this is the only file that needs to be modified (more instructions on how to modify it will follow);
- Resources: Directory with resources needed to run PIMBA, such as the adapters.txt file containing adapter sequences;
- Workflow: Directory with Snakemake scripts to run PIMBA;
- Test_data: Directory with data to test the algorithm;
- README.md: User guide similar to this document;
- pimba_smk_main.sh: Main script to run PIMBA after correctly editing the config.yaml file.

## How to Run PIMBA v3.0 in Snakemake?

### A) Configure the config.yaml File

The config.yaml file (inside the config folder) is the general configuration file. It should contain parameters such as the maximum number of processors, adapter sequences, input file paths, etc. Open the file in a text editor and make the following modifications. **Note: Use full paths; partial paths will not work in this version.**

#### General options for all modes

| Variable | Description |
| ----------- | ----------- |
| num_threads | The num_threads option indicates the maximum number of processors to be used for tasks that allow parallelization. |

#### General options for the Prepare Mode
If the user wishes to run the "prepare" mode to prepare the reads for the "run" mode, edit the following options:

| Variable | Description |
| ----------- | ----------- |
| minlength | The minimum length of a read after quality filtering. |
| minphred | The minimum PHRED score for quality filtering. |
| outputprepare | The name of the output file to be created in FASTA format. The .fasta extension is included automatically |

#### Inputs to run paired-end reads: 
If the user wishes to run PIMBA for paired-end reads, it is necessary to configure:

| Variable | Description |
| ----------- | ----------- |
| rawdatadir | the path to the directory where the reads are located |
| adapters | the path to the adapter file within the resources directory |

#### Inputs for single end and single index reads: 
If the user has unpaired reads containing a single index, the following options need to be configured:

| Variable | Description |
| ----------- | ----------- |
| raw_fastq | Input FASTQ file. |
| prefix | Name to be included as a prefix for the files to be generated. |
| barcodes_5end_txt | Path to the barcode file used as an index at the 5' end of the reads. |
| singleadapter | Adapter sequence found in the reads. |
| barcodes_5end_fasta | Path to the FASTA file for the barcodes_5end_txt. |

#### Inputs for single end and dual index reads: 
If the user has unpaired reads and a dual index, the following options need to be configured:

| Variable | Description |
| ----------- | ----------- |
|raw_fastq | Input FASTQ file. |
| barcodes_3end_txt | Path to the barcode file used as an index at the 3' end of the reads. |
| barcodes_3end_rev | Path to the reverse complement of the barcodes_3end_txt. |
| barcodes_3end_fasta | Path to the FASTA file for the barcodes_3end_txt. |
| barcodes_5end_dir | Path to the directory containing all the barcodes.fasta and barcodes.txt files used at the 5' end of the reads. Each 3' barcode should have a corresponding FASTA and TXT file with all the associated 5' barcodes. |
| forward_adapter | Forward primer sequence. |
| reverse_adapter | Reverse primer sequence. |
