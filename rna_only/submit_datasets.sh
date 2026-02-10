#!/bin/bash
# ============================================================================
# Multi-Dataset Pipeline Submission - RNA-only
# ============================================================================
#
# Submits the scMultiPreDICT RNA-only pipeline for MULTIPLE datasets.
#   - Within each dataset: steps run SERIALLY (SLURM dependency chains)
#   - Across datasets: pipelines run in PARALLEL (independent jobs)
#
# Setup:
#   1. Copy submit_settings.template.sh → submit_settings.sh, edit for cluster
#   2. Create per-dataset configs in datasets/ (copy from config_template.R)
#   3. Run: ./submit_datasets.sh
#
# Usage:
#   ./submit_datasets.sh                      # Submit ALL datasets in datasets/
#   ./submit_datasets.sh --dry-run            # Preview without submitting
#   ./submit_datasets.sh --datasets d1,d2     # Submit only specified datasets
#   ./submit_datasets.sh --start-step 2       # Start from step 2
#   ./submit_datasets.sh --stop-step 3        # Stop after step 3
#   ./submit_datasets.sh --help               # Show full usage
#
# Pipeline Steps (serial within each dataset):
#   1: Metacell Creation (PCA-based kNN smoothing)
#   2: Feature Extraction (HVG expression features)
#   3: Linear + Random Forest Model Training
#   4: Neural Network Training
#
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$SCRIPT_DIR"

# ─── Load cluster settings ──────────────────────────────────────────────────
SETTINGS_FILE="${PIPELINE_DIR}/submit_settings.sh"
if [ ! -f "$SETTINGS_FILE" ]; then
    echo "ERROR: Settings file not found: $SETTINGS_FILE"
    echo ""
    echo "Setup: cp submit_settings.template.sh submit_settings.sh"
    echo "Then edit submit_settings.sh with your cluster configuration."
    exit 1
fi
source "$SETTINGS_FILE"

# ─── Pipeline step definitions ───────────────────────────────────────────────
NUM_STEPS=4
declare -A STEP_NAMES=(
    [1]="metacell"
    [2]="feature_extract"
    [3]="linear_models"
    [4]="neural_network"
)
declare -A STEP_DESCRIPTIONS=(
    [1]="Metacell Creation (PCA + kNN smoothing)"
    [2]="Feature Extraction"
    [3]="Linear + Random Forest Models"
    [4]="Neural Network Training"
)

# Map step resources from settings
get_step_cpus()  { eval echo "\${STEP_${1}_CPUS}"; }
get_step_mem()   { eval echo "\${STEP_${1}_MEM}"; }
get_step_time()  { eval echo "\${STEP_${1}_TIME}"; }

# ─── Parse arguments ─────────────────────────────────────────────────────────
DRY_RUN=false
SPECIFIC_DATASETS=""
START_STEP=1
STOP_STEP=$NUM_STEPS

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true; shift ;;
        --datasets)
            SPECIFIC_DATASETS="$2"; shift 2 ;;
        --start-step)
            START_STEP="$2"; shift 2 ;;
        --stop-step)
            STOP_STEP="$2"; shift 2 ;;
        --help|-h)
            cat << 'EOF'
Multi-Dataset Pipeline Submission - RNA-only

Usage: ./submit_datasets.sh [OPTIONS]

Options:
  --dry-run             Preview jobs without submitting to SLURM
  --datasets d1,d2      Submit only named datasets (comma-separated, no .config.R)
  --start-step N        Start pipeline from step N (default: 1)
  --stop-step N         Stop pipeline after step N (default: 4)
  --help, -h            Show this help message

Steps:
  1  Metacell Creation      PCA-based kNN smoothing
  2  Feature Extraction     HVG expression features per target gene
  3  Linear Models          Ridge, Lasso, Elastic Net, Random Forest
  4  Neural Network         Feed-forward NN with grid search

Examples:
  # Submit all datasets, all steps
  ./submit_datasets.sh

  # Dry run to preview job chain
  ./submit_datasets.sh --dry-run

  # Submit only E7.5_rep1 and T_Cells starting from step 2
  ./submit_datasets.sh --datasets E7.5_rep1,T_Cells --start-step 2

  # Run only model training steps (3-4) for all datasets
  ./submit_datasets.sh --start-step 3 --stop-step 4
EOF
            exit 0 ;;
        *)
            echo "Unknown option: $1 (use --help for usage)"; exit 1 ;;
    esac
done

# ─── Discover dataset configs ────────────────────────────────────────────────
DATASET_DIR="${PIPELINE_DIR}/datasets"
if [ ! -d "$DATASET_DIR" ]; then
    echo "ERROR: Dataset directory not found: $DATASET_DIR"
    echo ""
    echo "Setup:"
    echo "  mkdir -p datasets"
    echo "  cp config_template.R datasets/MY_DATASET.config.R"
    echo "  # Edit datasets/MY_DATASET.config.R with dataset-specific paths"
    exit 1
fi

CONFIG_FILES=()
if [ -n "$SPECIFIC_DATASETS" ]; then
    IFS=',' read -ra DATASET_LIST <<< "$SPECIFIC_DATASETS"
    for ds in "${DATASET_LIST[@]}"; do
        cfg="${DATASET_DIR}/${ds}.config.R"
        if [ ! -f "$cfg" ]; then
            echo "ERROR: Config not found for dataset '${ds}': ${cfg}"
            exit 1
        fi
        CONFIG_FILES+=("$cfg")
    done
else
    shopt -s nullglob
    CONFIG_FILES=(${DATASET_DIR}/*.config.R)
    shopt -u nullglob
    if [ ${#CONFIG_FILES[@]} -eq 0 ]; then
        echo "ERROR: No .config.R files found in $DATASET_DIR"
        echo ""
        echo "Create dataset configs:"
        echo "  cp config_template.R datasets/MY_DATASET.config.R"
        exit 1
    fi
fi

# ─── Build active step list ──────────────────────────────────────────────────
ACTIVE_STEPS=()
for ((s=START_STEP; s<=STOP_STEP; s++)); do
    if [ "$s" -ge 1 ] && [ "$s" -le "$NUM_STEPS" ]; then
        ACTIVE_STEPS+=("$s")
    fi
done

if [ ${#ACTIVE_STEPS[@]} -eq 0 ]; then
    echo "ERROR: No valid steps in range $START_STEP-$STOP_STEP"
    exit 1
fi

# ─── Print submission plan ───────────────────────────────────────────────────
echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║       scMultiPreDICT - Multi-Dataset Submission             ║"
echo "║                   RNA-only Pipeline                         ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "  Pipeline dir:  $PIPELINE_DIR"
echo "  Partition:     $PARTITION"
echo "  Conda env:     $CONDA_ENV_R"
echo "  Log dir:       $LOG_BASE_DIR"
if $DRY_RUN; then
echo "  Mode:          DRY RUN (no jobs submitted)"
else
echo "  Mode:          LIVE SUBMISSION"
fi
echo ""
echo "  Datasets (${#CONFIG_FILES[@]}):"
for cfg in "${CONFIG_FILES[@]}"; do
    echo "    - $(basename "$cfg" .config.R)"
done
echo ""
echo "  Steps (${#ACTIVE_STEPS[@]}):"
for s in "${ACTIVE_STEPS[@]}"; do
    printf "    %d. %-35s [%s CPUs, %s RAM, %s]\n" \
        "$s" "${STEP_DESCRIPTIONS[$s]}" "$(get_step_cpus $s)" "$(get_step_mem $s)" "$(get_step_time $s)"
done
echo ""

if $DRY_RUN; then
    echo "─── DRY RUN: Showing commands that would be submitted ───"
    echo ""
fi

# ─── Submit job chains ───────────────────────────────────────────────────────
TOTAL_JOBS=0
declare -A ALL_JOB_CHAINS

for cfg in "${CONFIG_FILES[@]}"; do
    dataset_name=$(basename "$cfg" .config.R)
    config_path=$(realpath "$cfg")

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  Dataset: $dataset_name"
    echo "  Config:  $config_path"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    # Create log directory
    dataset_log_dir="${LOG_BASE_DIR}/${dataset_name}"
    mkdir -p "$dataset_log_dir"

    prev_job_id=""
    job_chain=""

    for step_num in "${ACTIVE_STEPS[@]}"; do
        step_name="${STEP_NAMES[$step_num]}"
        job_name="${step_name}_${dataset_name}"
        cpus=$(get_step_cpus $step_num)
        mem=$(get_step_mem $step_num)
        time_limit=$(get_step_time $step_num)

        # Build sbatch command
        sbatch_args=(
            --job-name="$job_name"
            -p "$PARTITION"
            -N 1
            --ntasks-per-node=1
            -c "$cpus"
            --mem="$mem"
            --time="$time_limit"
            --output="${dataset_log_dir}/${step_name}_%j.log"
            --error="${dataset_log_dir}/${step_name}_%j.err"
            --export="ALL,CONFIG_PATH=${config_path},STEP_NUM=${step_num},PIPELINE_DIR=${PIPELINE_DIR},CONDA_ENV_NAME=${CONDA_ENV_R}"
        )

        # Add dependency on previous step (serial within dataset)
        if [ -n "$prev_job_id" ]; then
            sbatch_args+=(--dependency="afterok:${prev_job_id}")
        fi

        sbatch_args+=("${PIPELINE_DIR}/slurm/run_step.sbatch")

        if $DRY_RUN; then
            echo "    [DRY RUN] sbatch ${sbatch_args[*]}"
            prev_job_id="DRY_${step_num}"
        else
            job_output=$(sbatch "${sbatch_args[@]}" 2>&1)
            if [ $? -ne 0 ]; then
                echo "    ERROR submitting step $step_num: $job_output"
                echo "    Aborting chain for dataset: $dataset_name"
                break
            fi
            prev_job_id=$(echo "$job_output" | grep -oP '\d+$')
            echo "    Step $step_num (${step_name}): Job ${prev_job_id}  [${cpus} CPUs, ${mem}, ${time_limit}]"
            TOTAL_JOBS=$((TOTAL_JOBS + 1))
        fi

        job_chain="${job_chain:+$job_chain → }${prev_job_id}"
    done

    echo "    Chain: $job_chain"
    echo ""
    ALL_JOB_CHAINS[$dataset_name]="$job_chain"
done

# ─── Print summary ───────────────────────────────────────────────────────────
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║                   SUBMISSION SUMMARY                        ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""

if $DRY_RUN; then
    echo "  DRY RUN complete. No jobs were submitted."
    echo "  Remove --dry-run to submit for real."
else
    echo "  Total jobs submitted: $TOTAL_JOBS"
    echo "  Datasets: ${#CONFIG_FILES[@]} (running in PARALLEL)"
    echo "  Steps per dataset: ${#ACTIVE_STEPS[@]} (running SERIALLY)"
    echo ""
    echo "  Job chains:"
    for dataset_name in "${!ALL_JOB_CHAINS[@]}"; do
        echo "    $dataset_name: ${ALL_JOB_CHAINS[$dataset_name]}"
    done
fi

echo ""
echo "  Monitor:  squeue -u $USER"
echo "  Cancel:   scancel <job_id>   (cancels downstream dependent jobs too)"
echo "  Logs:     $LOG_BASE_DIR/<dataset>/"
echo ""
