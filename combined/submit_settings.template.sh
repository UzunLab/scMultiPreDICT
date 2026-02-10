#!/bin/bash
# ============================================================================
# Cluster Settings Template - Combined (RNA+ATAC) Pipeline
# ============================================================================
# Copy to submit_settings.sh and edit for your cluster:
#   cp submit_settings.template.sh submit_settings.sh
# ============================================================================

# ─── SLURM Partition ─────────────────────────────────────────────────────────
PARTITION="compute"
GPU_PARTITION="gpu"

# ─── Conda Environment Names ─────────────────────────────────────────────────
CONDA_ENV_R="r-bioc-43"
CONDA_ENV_PYTHON="python_env"

# ─── Log Directory ───────────────────────────────────────────────────────────
# Logs are organized as: LOG_BASE_DIR/<dataset_name>/<step_name>_<jobid>.log
LOG_BASE_DIR="$HOME/scMultiPreDICT_logs/combined"

# ─── Per-Step Resource Allocation ─────────────────────────────────────────────

# Step 1: Quality Control & Preprocessing
STEP_1_CPUS=8
STEP_1_MEM="64G"
STEP_1_TIME="04:00:00"

# Step 2: Data Splitting
STEP_2_CPUS=4
STEP_2_MEM="32G"
STEP_2_TIME="01:00:00"

# Step 3: Metacell Creation (method-dependent)
STEP_3_CPUS=32
STEP_3_MEM="256G"
STEP_3_TIME="08:00:00"

# Step 4: Feature Extraction
STEP_4_CPUS=32
STEP_4_MEM="256G"
STEP_4_TIME="12:00:00"

# Step 5: Linear + Random Forest Model Training
STEP_5_CPUS=32
STEP_5_MEM="128G"
STEP_5_TIME="12:00:00"

# Step 6: Neural Network Training
STEP_6_CPUS=16
STEP_6_MEM="64G"
STEP_6_TIME="24:00:00"

# ─── Optional: Python Autoencoder Step (for multivi/scvi_peakvi methods) ────
# Only used when --include-python is passed to submit_datasets.sh
PYTHON_STEP_CPUS=8
PYTHON_STEP_MEM="64G"
PYTHON_STEP_TIME="06:00:00"
PYTHON_STEP_GPUS=1
