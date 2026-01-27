#!/bin/bash
#SBATCH --job-name=model_nn_YOUR_SAMPLE_NAME
#SBATCH -p compute
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/Home/evm5863/scATAC_only/jobs/logs/model_nn/%x_%j.log
#SBATCH --error=/gpfs/Home/evm5863/scATAC_only/jobs/logs/model_nn/%x_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH --mem=200G

set -euo pipefail

# Initialize and activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate r-bioc-43

# Set R library paths (conda R)
export R_LIBS="$CONDA_PREFIX/lib/R/library"
export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
export R_LIBS_SITE=""

# Run the script
Rscript 06_neural_network.R 

echo "Job completed at $(date)"
echo "Exit code: $?"
