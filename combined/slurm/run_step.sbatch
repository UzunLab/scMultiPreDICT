#!/bin/bash
# ============================================================================
# Generic Pipeline Step Runner - Combined (RNA+ATAC)
# ============================================================================
# This template is submitted by submit_datasets.sh with environment variables.
# Do NOT submit this directly - use submit_datasets.sh instead.
#
# Required environment variables (set via sbatch --export):
#   CONFIG_PATH    - Absolute path to the dataset's .config.R file
#   STEP_NUM       - Pipeline step number to run (1-6)
#   PIPELINE_DIR   - Absolute path to the combined/ directory
#   CONDA_ENV_NAME - Conda environment name for R
#
# Pipeline Steps:
#   1: Quality Control
#   2: Data Splitting
#   3: Metacell Creation (method depends on config)
#   4: Feature Extraction
#   5: Linear and Tree-Based Model Training
#   6: Neural Network Training
# ============================================================================

set -euo pipefail

# ─── Validate required environment variables ─────────────────────────────────
for var in CONFIG_PATH STEP_NUM PIPELINE_DIR CONDA_ENV_NAME; do
    if [ -z "${!var:-}" ]; then
        echo "ERROR: Required environment variable $var is not set."
        echo "This script should be submitted via submit_datasets.sh"
        exit 1
    fi
done

# ─── Initialize conda ────────────────────────────────────────────────────────
CONDA_FOUND=false
for conda_init in "$HOME/miniconda3/etc/profile.d/conda.sh" \
                  "$HOME/anaconda3/etc/profile.d/conda.sh" \
                  "$HOME/miniforge3/etc/profile.d/conda.sh" \
                  "$HOME/mambaforge/etc/profile.d/conda.sh"; do
    if [ -f "$conda_init" ]; then
        source "$conda_init"
        CONDA_FOUND=true
        break
    fi
done

if ! $CONDA_FOUND; then
    echo "ERROR: Cannot find conda initialization script"
    echo "Searched: ~/miniconda3, ~/anaconda3, ~/miniforge3, ~/mambaforge"
    exit 1
fi

conda activate "$CONDA_ENV_NAME"

# Set R library paths for conda R
export R_LIBS="$CONDA_PREFIX/lib/R/library"
export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
export R_LIBS_SITE=""

# ─── Print job information ───────────────────────────────────────────────────
echo "============================================"
echo "scMultiPreDICT - Combined Pipeline Step"
echo "============================================"
echo "Job ID:        ${SLURM_JOB_ID:-local}"
echo "Job Name:      ${SLURM_JOB_NAME:-local}"
echo "Config:        ${CONFIG_PATH}"
echo "Step:          ${STEP_NUM}"
echo "Pipeline Dir:  ${PIPELINE_DIR}"
echo "Conda Env:     ${CONDA_DEFAULT_ENV:-$CONDA_ENV_NAME}"
echo "CPUs:          ${SLURM_CPUS_PER_TASK:-N/A}"
echo "Memory:        ${SLURM_MEM_PER_NODE:-N/A}MB"
echo "Node:          ${SLURM_NODELIST:-local}"
echo "Start Time:    $(date)"
echo "============================================"

# ─── Run the pipeline step ───────────────────────────────────────────────────
cd "$PIPELINE_DIR"

echo ""
echo "Running: Rscript run_pipeline.R --config \"$CONFIG_PATH\" --steps $STEP_NUM"
echo ""

Rscript run_pipeline.R --config "$CONFIG_PATH" --steps "$STEP_NUM"

EXIT_CODE=$?

echo ""
echo "============================================"
echo "Step $STEP_NUM completed at $(date)"
echo "Exit code: $EXIT_CODE"
echo "============================================"

exit $EXIT_CODE
