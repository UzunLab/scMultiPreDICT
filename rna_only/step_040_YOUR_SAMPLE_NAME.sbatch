#!/bin/bash
#SBATCH --job-name=feature_extract_YOUR_SAMPLE_NAME
#SBATCH -p compute
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/Home/evm5863/scRNA_only/jobs/logs/feature_extract/%x_%j.log
#SBATCH --error=/gpfs/Home/evm5863/scRNA_only/jobs/logs/feature_extract/%x_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 32
#SBATCH --mem=128G

set -euo pipefail

# Create log directory
mkdir -p "/gpfs/Home/evm5863/scRNA_only/jobs/logs/feature_extract"

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
Rscript 04_feature_extraction.R 

# Print completion
echo ""
echo "Job completed at $(date)"
echo "Exit code: $?"
