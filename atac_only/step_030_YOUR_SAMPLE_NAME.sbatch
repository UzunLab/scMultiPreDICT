#!/bin/bash
#SBATCH --job-name=metacell_YOUR_SAMPLE_NAME
#SBATCH -p compute
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/Home/evm5863/scATAC_only/jobs/logs/metacell/%x_%j.log
#SBATCH --error=/gpfs/Home/evm5863/scATAC_only/jobs/logs/metacell/%x_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 32
#SBATCH --mem=128G

set -euo pipefail

# Initialize and activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate r-bioc-43

# Set R library paths (conda R)
export R_LIBS="$CONDA_PREFIX/lib/R/library"
export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
export R_LIBS_SITE=""

# Run the script
Rscript 03_metacell_creation.R 

echo "Job completed at $(date)"
echo "Exit code: $?"
