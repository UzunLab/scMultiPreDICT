#!/bin/bash
# ============================================================================
# Python Step Runner - Combined Pipeline
# ============================================================================
# For running Python-based steps (e.g., autoencoder training for multivi/scvi_peakvi).
# Submitted by submit_datasets.sh when INCLUDE_PYTHON_STEPS is enabled.
#
# Required environment variables (set via sbatch --export):
#   CONFIG_PATH    - Absolute path to the dataset's .config.R file
#   PYTHON_CMD     - Full Python command to execute
#   PIPELINE_DIR   - Absolute path to the combined/ directory
#   CONDA_ENV_NAME - Conda environment name for Python
# ============================================================================

set -euo pipefail

# ─── Validate required environment variables ─────────────────────────────────
for var in PYTHON_CMD PIPELINE_DIR CONDA_ENV_NAME; do
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
    exit 1
fi

conda activate "$CONDA_ENV_NAME"

# ─── Print job information ───────────────────────────────────────────────────
echo "============================================"
echo "scMultiPreDICT - Python Step"
echo "============================================"
echo "Job ID:        ${SLURM_JOB_ID:-local}"
echo "Job Name:      ${SLURM_JOB_NAME:-local}"
echo "Config:        ${CONFIG_PATH:-N/A}"
echo "Command:       ${PYTHON_CMD}"
echo "Pipeline Dir:  ${PIPELINE_DIR}"
echo "Conda Env:     ${CONDA_DEFAULT_ENV:-$CONDA_ENV_NAME}"
echo "Start Time:    $(date)"
echo "============================================"

# ─── Run the Python command ──────────────────────────────────────────────────
cd "$PIPELINE_DIR"

echo ""
echo "Running: $PYTHON_CMD"
echo ""

eval "$PYTHON_CMD"

EXIT_CODE=$?

echo ""
echo "============================================"
echo "Python step completed at $(date)"
echo "Exit code: $EXIT_CODE"
echo "============================================"

exit $EXIT_CODE
