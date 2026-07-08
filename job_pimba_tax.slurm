#!/bin/bash
#SBATCH --job-name=pb_tax_16s_rdp      # Job name (customize as needed; e.g., pb_tax_16s_silva, pb_tax_coi_bold_ncbi)
#SBATCH --output=%x_%j.out               # Standard output file
#SBATCH --error=%x_%j.err                # Standard error file
#SBATCH --partition=op1_slurm            # Partition or queue name
#SBATCH --time=03:00:00                  # Max time
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks=1                       # Number of tasks per node
#SBATCH --cpus-per-task=8                # Number of cpus per task
#SBATCH --mem=20G                        # Max RAM usage

printf "Begin: `date`\n";
TBEGIN=`echo "print time();" | perl`

printf "\n"

cd $SLURM_SUBMIT_DIR

# Load snakemake module
module load snakemake

srun bash pimba_smk_tax.sh \
    -r 16S-RDP \
    -t "$SLURM_CPUS_PER_TASK" \
    -c config/config_tax.yaml \
    -d .

TEND=`echo "print time();" | perl`

printf "\n"
printf "End: `date`\n";
printf "Run time (s): `expr $TEND - $TBEGIN`\n";
printf "Run time (min): `expr $(( ($TEND - $TBEGIN)/60 ))`\n";