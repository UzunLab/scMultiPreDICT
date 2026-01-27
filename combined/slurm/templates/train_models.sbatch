#!/bin/bash
#SBATCH --job-name=model_linear_E7.5_rep2
#SBATCH -p compute
#SBATCH --time=12:00:00
#SBATCH --output=/logs/model_linear/%x_%j.log
#SBATCH --error=/logs/model_linear/%x_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 32
#SBATCH --mem=128G

set -euo pipefail

# Create log directory
mkdir -p "/logs/model_linear"

# Initialize and activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate r-bioc-43

# Print environment info
echo "============================================"
echo "Job Information"
echo "============================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "CPUs: 32"
echo "Memory: 128G"
echo "Conda Env: $CONDA_DEFAULT_ENV"
echo "Start Time: $(date)"
echo "============================================"

# Set R library paths
export R_LIBS="$CONDA_PREFIX/lib/R/library"
export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
export R_LIBS_SITE=""

# Run the script
Rscript 05_linear_tree_models.R 

# Print completion
echo ""
echo "Job completed at $(date)"
echo "Exit code: $?"
