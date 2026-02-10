#!/bin/bash
# ============================================================================
# Cluster Settings Template - ATAC-only Pipeline
# ============================================================================
# Copy to submit_settings.sh and edit for your cluster:
#   cp submit_settings.template.sh submit_settings.sh
# ============================================================================

# ─── SLURM Partition ─────────────────────────────────────────────────────────
PARTITION="compute"

# ─── Conda Environment Names ─────────────────────────────────────────────────
CONDA_ENV_R="r-bioc-43"

# ─── Log Directory ───────────────────────────────────────────────────────────
# Logs are organized as: LOG_BASE_DIR/<dataset_name>/<step_name>_<jobid>.log
LOG_BASE_DIR="$HOME/scMultiPreDICT_logs/atac_only"

# ─── Per-Step Resource Allocation ─────────────────────────────────────────────
# Step 1: Metacell Creation (LSI-based kNN smoothing)
STEP_1_CPUS=32
STEP_1_MEM="128G"
STEP_1_TIME="06:00:00"

# Step 2: Feature Extraction (peak accessibility features)
STEP_2_CPUS=32
STEP_2_MEM="128G"
STEP_2_TIME="08:00:00"

# Step 3: Linear + Random Forest Model Training
STEP_3_CPUS=32
STEP_3_MEM="200G"
STEP_3_TIME="12:00:00"

# Step 4: Neural Network Training
STEP_4_CPUS=16
STEP_4_MEM="200G"
STEP_4_TIME="24:00:00"
