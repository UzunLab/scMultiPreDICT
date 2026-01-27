#!/bin/bash
#SBATCH --job-name=model_nn_YOUR_SAMPLE_NAME
#SBATCH -p compute
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/Home/evm5863/scRNA_only/jobs/logs/model_nn/%x_%j.log
#SBATCH --error=/gpfs/Home/evm5863/scRNA_only/jobs/logs/model_nn/%x_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH --mem=200G

set -euo pipefail

# Create log directory
mkdir -p "/gpfs/Home/evm5863/scRNA_only/jobs/logs/model_nn"

# Initialize and activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate r-bioc-43

# Print environment info
echo "============================================"
echo "Job Information"
echo "============================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "CPUs: 16"
echo "Memory: 200G"
echo "Conda Env: $CONDA_DEFAULT_ENV"
echo "Start Time: $(date)"
echo "============================================"

# Set R library paths
export R_LIBS="$CONDA_PREFIX/lib/R/library"
export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
export R_LIBS_SITE=""

# Run the script
Rscript 06_neural_network.R 

# Print completion
echo ""
echo "Job completed at $(date)"
echo "Exit code: $?"
