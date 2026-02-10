#!/bin/bash
# ============================================================================
# Master Orchestration - All Three Pipelines (Combined → RNA + ATAC)
# ============================================================================
#
# Submits the full scMultiPreDICT analysis across all three pipelines with
# proper cross-pipeline dependency handling:
#
#   COMBINED:  step1 → step2 → step3 → step4 → step5 → step6
#                         ↓ (after step 2 completes)
#   RNA_ONLY:           rna1 → rna2 → rna3 → rna4
#                         ↓ (after step 2 completes)
#   ATAC_ONLY:         atac1 → atac2 → atac3 → atac4
#
# Combined step 2 produces the Seurat splits that
# RNA-only and ATAC-only pipelines consume. After step 2 finishes,
# three independent branches run in parallel:
#   - Combined steps 3-6 (combined metacell → NN)
#   - RNA-only steps 1-4  (RNA metacell → NN)
#   - ATAC-only steps 1-4 (ATAC metacell → NN)
#
# Setup (one-time):
#   1. Configure each pipeline's cluster settings:
#        cp combined/submit_settings.template.sh combined/submit_settings.sh
#        cp rna_only/submit_settings.template.sh  rna_only/submit_settings.sh
#        cp atac_only/submit_settings.template.sh atac_only/submit_settings.sh
#   2. Create per-dataset configs for the combined pipeline:
#        cp combined/R/config_template.R combined/datasets/T_Cells.config.R
#   3. Create matching configs for RNA/ATAC (optional — only needed pipelines):
#        cp rna_only/config_template.R  rna_only/datasets/T_Cells.config.R
#        cp atac_only/config_template.R atac_only/datasets/T_Cells.config.R
#   4. Edit all config files with dataset-specific paths
#   5. Run: ./submit_all.sh
#
# Usage:
#   ./submit_all.sh                           # All datasets, all pipelines
#   ./submit_all.sh --dry-run                 # Preview without submitting
#   ./submit_all.sh --datasets d1,d2          # Specific datasets only
#   ./submit_all.sh --pipelines combined,rna  # Skip ATAC pipeline
#   ./submit_all.sh --help                    # Full usage
#
# ============================================================================

set -eo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ─── Defaults ─────────────────────────────────────────────────────────────────
DRY_RUN=false
SPECIFIC_DATASETS=""
REQUESTED_PIPELINES="combined,rna,atac"  # all by default
DRY_COUNTER_FILE=$(mktemp)
echo "0" > "$DRY_COUNTER_FILE"
trap "rm -f '$DRY_COUNTER_FILE'" EXIT

# ─── Parse arguments ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true; shift ;;
        --datasets)
            SPECIFIC_DATASETS="$2"; shift 2 ;;
        --pipelines)
            REQUESTED_PIPELINES="$2"; shift 2 ;;
        --help|-h)
            cat << 'EOF'
Master Orchestration - scMultiPreDICT All Pipelines

Usage: ./submit_all.sh [OPTIONS]

Options:
  --dry-run              Preview jobs without submitting to SLURM
  --datasets d1,d2       Submit only named datasets (comma-separated)
  --pipelines p1,p2      Which pipelines to run (default: combined,rna,atac)
  --help, -h             Show this help message

Pipelines (use short names with --pipelines):
  combined   Combined RNA+ATAC pipeline (6 steps) — always runs first
  rna        RNA-only pipeline (4 steps) — depends on combined step 2
  atac       ATAC-only pipeline (4 steps) — depends on combined step 2

Dependency Graph (per dataset):
  COMBINED:  step1 → step2 → step3 → step4 → step5 → step6
                        ↓ (afterok)
  RNA_ONLY:           rna1 → rna2 → rna3 → rna4
                        ↓ (afterok)
  ATAC_ONLY:         atac1 → atac2 → atac3 → atac4

Dataset configs:
  Each pipeline reads from <pipeline>/datasets/<name>.config.R
  The combined pipeline config is REQUIRED for each dataset.
  RNA-only and ATAC-only configs are optional — if missing, that
  pipeline is skipped for that dataset.

Examples:
  # Full pipeline for all datasets
  ./submit_all.sh

  # Preview everything
  ./submit_all.sh --dry-run

  # Only T_Cells and E7.5, combined + RNA only
  ./submit_all.sh --datasets T_Cells,E7.5_rep1 --pipelines combined,rna

  # Only combined pipeline (same as cd combined && ./submit_datasets.sh)
  ./submit_all.sh --pipelines combined
EOF
            exit 0 ;;
        *)
            echo "Unknown option: $1 (use --help for usage)"; exit 1 ;;
    esac
done

# ─── Parse requested pipelines ───────────────────────────────────────────────
RUN_COMBINED=false
RUN_RNA=false
RUN_ATAC=false

IFS=',' read -ra PIPE_LIST <<< "$REQUESTED_PIPELINES"
for p in "${PIPE_LIST[@]}"; do
    case "$p" in
        combined) RUN_COMBINED=true ;;
        rna)      RUN_RNA=true ;;
        atac)     RUN_ATAC=true ;;
        *)        echo "ERROR: Unknown pipeline '$p'. Use: combined, rna, atac"; exit 1 ;;
    esac
done

if ! $RUN_COMBINED && ($RUN_RNA || $RUN_ATAC); then
    echo "WARNING: RNA-only and ATAC-only pipelines depend on the combined pipeline."
    echo "         Running without combined assumes combined steps 1-2 are already done."
    echo ""
fi

# ─── Load cluster settings from each pipeline ────────────────────────────────
# We namespace all variables to avoid collisions between pipelines.

load_settings() {
    local p="$1" pfx="$2"
    local f="${REPO_DIR}/${p}/submit_settings.sh"
    if [ ! -f "$f" ]; then
        echo "ERROR: Settings file not found: $f"
        echo "  Setup: cp ${p}/submit_settings.template.sh ${p}/submit_settings.sh"
        exit 1
    fi
    # Clear previous values before sourcing
    unset PARTITION CONDA_ENV_R LOG_BASE_DIR 2>/dev/null || true
    for i in 1 2 3 4 5 6; do
        unset "STEP_${i}_CPUS" "STEP_${i}_MEM" "STEP_${i}_TIME" 2>/dev/null || true
    done
    source "$f"
    # Copy into prefixed variables
    eval "${pfx}_PARTITION=\"$PARTITION\""
    eval "${pfx}_CONDA_ENV_R=\"$CONDA_ENV_R\""
    eval "${pfx}_LOG_BASE_DIR=\"$LOG_BASE_DIR\""
    for i in 1 2 3 4 5 6; do
        local cv="STEP_${i}_CPUS" mv="STEP_${i}_MEM" tv="STEP_${i}_TIME"
        if [ -n "${!cv:-}" ]; then
            eval "${pfx}_STEP_${i}_CPUS=\"${!cv}\""
            eval "${pfx}_STEP_${i}_MEM=\"${!mv}\""
            eval "${pfx}_STEP_${i}_TIME=\"${!tv}\""
        fi
    done
}

if $RUN_COMBINED; then load_settings "combined" "COMB"; fi
if $RUN_RNA;      then load_settings "rna_only" "RNA"; fi
if $RUN_ATAC;     then load_settings "atac_only" "ATAC"; fi

# ─── Pipeline step definitions ───────────────────────────────────────────────
# Combined pipeline: 6 steps
declare -A COMB_STEP_NAMES=(
    [1]="qc_preprocess" [2]="data_split" [3]="metacell"
    [4]="feature_extract" [5]="linear_models" [6]="neural_network"
)
declare -A COMB_STEP_DESC=(
    [1]="Quality Control & Preprocessing"
    [2]="Data Splitting + Target Gene Selection"
    [3]="Metacell Creation (method from config)"
    [4]="Feature Extraction"
    [5]="Linear + RF Model Training & Prediction"
    [6]="Neural Network Training & Prediction"
)

# RNA-only pipeline: 4 steps
declare -A RNA_STEP_NAMES=(
    [1]="metacell" [2]="feature_extract" [3]="linear_models" [4]="neural_network"
)
declare -A RNA_STEP_DESC=(
    [1]="Metacell Creation (PCA + kNN)"
    [2]="Feature Extraction"
    [3]="Linear + RF Model Training & Prediction"
    [4]="Neural Network Training & Prediction"
)

# ATAC-only pipeline: 4 steps
declare -A ATAC_STEP_NAMES=(
    [1]="metacell" [2]="feature_extract" [3]="linear_models" [4]="neural_network"
)
declare -A ATAC_STEP_DESC=(
    [1]="Metacell Creation (LSI + kNN)"
    [2]="Feature Extraction (peak accessibility)"
    [3]="Linear + RF Model Training & Prediction"
    [4]="Neural Network Training & Prediction"
)

# ─── Helper: get step resources ──────────────────────────────────────────────
get_resource() {
    local prefix="$1" step="$2" resource="$3"
    eval echo "\${${prefix}_STEP_${step}_${resource}}"
}

# ─── Helper: submit a single SLURM step ─────────────────────────────────────
# Sets RETURNED_JOB_ID to the job ID (or DRY_N in dry-run mode)
RETURNED_JOB_ID=""

submit_step() {
    local job_name="$1"
    local partition="$2"
    local cpus="$3"
    local mem="$4"
    local time_limit="$5"
    local log_dir="$6"
    local log_prefix="$7"
    local conda_env="$8"
    local config_path="$9"
    local step_num="${10}"
    local pipeline_dir="${11}"
    local sbatch_template="${12}"
    local prev_job_id="${13:-}"

    local sbatch_args=(
        --job-name="$job_name"
        -p "$partition"
        -N 1 --ntasks-per-node=1
        -c "$cpus"
        --mem="$mem"
        --time="$time_limit"
        --output="${log_dir}/${log_prefix}_%j.log"
        --error="${log_dir}/${log_prefix}_%j.err"
        --export="ALL,CONFIG_PATH=${config_path},STEP_NUM=${step_num},PIPELINE_DIR=${pipeline_dir},CONDA_ENV_NAME=${conda_env}"
    )

    if [ -n "$prev_job_id" ]; then
        sbatch_args+=(--dependency="afterok:${prev_job_id}")
    fi

    sbatch_args+=("$sbatch_template")

    if $DRY_RUN; then
        local cnt
        cnt=$(cat "$DRY_COUNTER_FILE")
        cnt=$((cnt + 1))
        echo "$cnt" > "$DRY_COUNTER_FILE"
        echo "    [DRY RUN] sbatch ${sbatch_args[*]}"
        RETURNED_JOB_ID="DRY_${cnt}"
    else
        local job_output
        job_output=$(sbatch "${sbatch_args[@]}" 2>&1)
        if [ $? -ne 0 ]; then
            echo "    ERROR: $job_output"
            RETURNED_JOB_ID="FAILED"
            return 1
        fi
        RETURNED_JOB_ID=$(echo "$job_output" | grep -oP '\d+$')
        echo "    Submitted: Job ${RETURNED_JOB_ID}  [${cpus} CPUs, ${mem}, ${time_limit}]"
    fi
}

# ─── Discover datasets ──────────────────────────────────────────────────────
# Datasets are discovered from combined/datasets/ (required)
COMBINED_DATASET_DIR="${REPO_DIR}/combined/datasets"
RNA_DATASET_DIR="${REPO_DIR}/rna_only/datasets"
ATAC_DATASET_DIR="${REPO_DIR}/atac_only/datasets"

if $RUN_COMBINED && [ ! -d "$COMBINED_DATASET_DIR" ]; then
    echo "ERROR: Combined dataset directory not found: $COMBINED_DATASET_DIR"
    exit 1
fi

# Build dataset list
DATASETS=()
if [ -n "$SPECIFIC_DATASETS" ]; then
    IFS=',' read -ra DATASETS <<< "$SPECIFIC_DATASETS"
else
    # Auto-discover from combined/datasets/
    if $RUN_COMBINED; then
        shopt -s nullglob
        for cfg_file in "${COMBINED_DATASET_DIR}"/*.config.R; do
            DATASETS+=("$(basename "$cfg_file" .config.R)")
        done
        shopt -u nullglob
    else
        # If not running combined, discover from whichever pipeline is requested
        if $RUN_RNA; then
            shopt -s nullglob
            for cfg_file in "${RNA_DATASET_DIR}"/*.config.R; do
                DATASETS+=("$(basename "$cfg_file" .config.R)")
            done
            shopt -u nullglob
        elif $RUN_ATAC; then
            shopt -s nullglob
            for cfg_file in "${ATAC_DATASET_DIR}"/*.config.R; do
                DATASETS+=("$(basename "$cfg_file" .config.R)")
            done
            shopt -u nullglob
        fi
    fi
fi

if [ ${#DATASETS[@]} -eq 0 ]; then
    echo "ERROR: No datasets found."
    if $RUN_COMBINED; then
        echo "  Create configs in combined/datasets/:"
        echo "    cp combined/R/config_template.R combined/datasets/MY_DATASET.config.R"
    fi
    exit 1
fi

# ─── Validate configs exist ─────────────────────────────────────────────────
echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║       scMultiPreDICT - Master Pipeline Orchestration        ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
if $DRY_RUN; then
    echo "  Mode:     DRY RUN (no jobs submitted)"
else
    echo "  Mode:     LIVE SUBMISSION"
fi
echo "  Repo:     $REPO_DIR"
echo ""

# For each dataset, determine which pipelines have configs
declare -A DS_HAS_COMBINED
declare -A DS_HAS_RNA
declare -A DS_HAS_ATAC

echo "  Datasets and available pipelines:"
for ds in "${DATASETS[@]}"; do
    pipelines_available=""

    if $RUN_COMBINED; then
        cfg="${COMBINED_DATASET_DIR}/${ds}.config.R"
        if [ -f "$cfg" ]; then
            DS_HAS_COMBINED[$ds]="$cfg"
            pipelines_available="${pipelines_available} combined"
        else
            echo "  ERROR: Combined config missing for '$ds': $cfg"
            exit 1
        fi
    fi

    if $RUN_RNA; then
        cfg="${RNA_DATASET_DIR}/${ds}.config.R"
        if [ -f "$cfg" ]; then
            DS_HAS_RNA[$ds]="$cfg"
            pipelines_available="${pipelines_available} rna"
        else
            echo "    NOTE: No RNA-only config for '$ds' — skipping RNA pipeline"
        fi
    fi

    if $RUN_ATAC; then
        cfg="${ATAC_DATASET_DIR}/${ds}.config.R"
        if [ -f "$cfg" ]; then
            DS_HAS_ATAC[$ds]="$cfg"
            pipelines_available="${pipelines_available} atac"
        else
            echo "    NOTE: No ATAC-only config for '$ds' — skipping ATAC pipeline"
        fi
    fi

    echo "    ${ds}: [${pipelines_available# } ]"
done

echo ""

# ─── Show step plan per pipeline ─────────────────────────────────────────────
if $RUN_COMBINED; then
    echo "  Combined pipeline (6 steps):"
    echo "    Partition: $COMB_PARTITION | Conda: $COMB_CONDA_ENV_R | Logs: $COMB_LOG_BASE_DIR"
    for s in 1 2 3 4 5 6; do
        printf "    %d. %-42s [%s CPUs, %s RAM, %s]\n" \
            "$s" "${COMB_STEP_DESC[$s]}" \
            "$(get_resource COMB $s CPUS)" "$(get_resource COMB $s MEM)" "$(get_resource COMB $s TIME)"
    done
    echo ""
fi

if $RUN_RNA && [ ${#DS_HAS_RNA[@]} -gt 0 ]; then
    echo "  RNA-only pipeline (4 steps, starts after combined step 2):"
    echo "    Partition: $RNA_PARTITION | Conda: $RNA_CONDA_ENV_R | Logs: $RNA_LOG_BASE_DIR"
    for s in 1 2 3 4; do
        printf "    %d. %-42s [%s CPUs, %s RAM, %s]\n" \
            "$s" "${RNA_STEP_DESC[$s]}" \
            "$(get_resource RNA $s CPUS)" "$(get_resource RNA $s MEM)" "$(get_resource RNA $s TIME)"
    done
    echo ""
fi

if $RUN_ATAC && [ ${#DS_HAS_ATAC[@]} -gt 0 ]; then
    echo "  ATAC-only pipeline (4 steps, starts after combined step 2):"
    echo "    Partition: $ATAC_PARTITION | Conda: $ATAC_CONDA_ENV_R | Logs: $ATAC_LOG_BASE_DIR"
    for s in 1 2 3 4; do
        printf "    %d. %-42s [%s CPUs, %s RAM, %s]\n" \
            "$s" "${ATAC_STEP_DESC[$s]}" \
            "$(get_resource ATAC $s CPUS)" "$(get_resource ATAC $s MEM)" "$(get_resource ATAC $s TIME)"
    done
    echo ""
fi

if $DRY_RUN; then
    echo "─── DRY RUN: Showing commands that would be submitted ───"
    echo ""
fi

# ─── Submit job chains per dataset ───────────────────────────────────────────
TOTAL_JOBS=0

for ds in "${DATASETS[@]}"; do
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  Dataset: $ds"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    combined_step2_job=""
    combined_chain=""
    rna_chain=""
    atac_chain=""

    # ── Combined pipeline ─────────────────────────────────────────────────
    if $RUN_COMBINED && [ -n "${DS_HAS_COMBINED[$ds]:-}" ]; then
        config_path=$(realpath "${DS_HAS_COMBINED[$ds]}")
        log_dir="${COMB_LOG_BASE_DIR}/${ds}"
        mkdir -p "$log_dir"

        echo "  ┌─ Combined Pipeline"

        prev_job=""
        for step in 1 2 3 4 5 6; do
            step_name="${COMB_STEP_NAMES[$step]}"
            job_name="comb_${step_name}_${ds}"
            cpus=$(get_resource COMB $step CPUS)
            mem=$(get_resource COMB $step MEM)
            time_limit=$(get_resource COMB $step TIME)

            echo "  │  Step $step: ${COMB_STEP_DESC[$step]}"

            submit_step \
                "$job_name" "$COMB_PARTITION" "$cpus" "$mem" "$time_limit" \
                "$log_dir" "$step_name" "$COMB_CONDA_ENV_R" \
                "$config_path" "$step" "${REPO_DIR}/combined" \
                "${REPO_DIR}/combined/slurm/run_step.sbatch" \
                "$prev_job"

            if [ "$RETURNED_JOB_ID" = "FAILED" ]; then
                echo "  │  ✗ Aborting combined chain for $ds"
                break
            fi

            combined_chain="${combined_chain:+$combined_chain → }${RETURNED_JOB_ID}"
            prev_job="$RETURNED_JOB_ID"
            TOTAL_JOBS=$((TOTAL_JOBS + 1))

            # Capture step 2 job ID for rna/atac dependency
            if [ "$step" -eq 2 ]; then
                combined_step2_job="$RETURNED_JOB_ID"
            fi
        done

        echo "  │  Chain: $combined_chain"
        echo "  └──────────────────────────"
    fi

    # ── RNA-only pipeline (depends on combined step 2) ────────────────────
    if $RUN_RNA && [ -n "${DS_HAS_RNA[$ds]:-}" ]; then
        config_path=$(realpath "${DS_HAS_RNA[$ds]}")
        log_dir="${RNA_LOG_BASE_DIR}/${ds}"
        mkdir -p "$log_dir"

        echo "  ┌─ RNA-only Pipeline"
        if [ -n "$combined_step2_job" ]; then
            echo "  │  (depends on combined step 2: job $combined_step2_job)"
        fi

        prev_job="$combined_step2_job"  # first RNA step depends on combined step 2
        for step in 1 2 3 4; do
            step_name="${RNA_STEP_NAMES[$step]}"
            job_name="rna_${step_name}_${ds}"
            cpus=$(get_resource RNA $step CPUS)
            mem=$(get_resource RNA $step MEM)
            time_limit=$(get_resource RNA $step TIME)

            echo "  │  Step $step: ${RNA_STEP_DESC[$step]}"

            submit_step \
                "$job_name" "$RNA_PARTITION" "$cpus" "$mem" "$time_limit" \
                "$log_dir" "$step_name" "$RNA_CONDA_ENV_R" \
                "$config_path" "$step" "${REPO_DIR}/rna_only" \
                "${REPO_DIR}/rna_only/slurm/run_step.sbatch" \
                "$prev_job"

            if [ "$RETURNED_JOB_ID" = "FAILED" ]; then
                echo "  │  ✗ Aborting RNA chain for $ds"
                break
            fi

            rna_chain="${rna_chain:+$rna_chain → }${RETURNED_JOB_ID}"
            prev_job="$RETURNED_JOB_ID"
            TOTAL_JOBS=$((TOTAL_JOBS + 1))
        done

        echo "  │  Chain: ${combined_step2_job:+${combined_step2_job} (combined step 2) → }$rna_chain"
        echo "  └──────────────────────────"
    fi

    # ── ATAC-only pipeline (depends on combined step 2) ───────────────────
    if $RUN_ATAC && [ -n "${DS_HAS_ATAC[$ds]:-}" ]; then
        config_path=$(realpath "${DS_HAS_ATAC[$ds]}")
        log_dir="${ATAC_LOG_BASE_DIR}/${ds}"
        mkdir -p "$log_dir"

        echo "  ┌─ ATAC-only Pipeline"
        if [ -n "$combined_step2_job" ]; then
            echo "  │  (depends on combined step 2: job $combined_step2_job)"
        fi

        prev_job="$combined_step2_job"  # first ATAC step depends on combined step 2
        for step in 1 2 3 4; do
            step_name="${ATAC_STEP_NAMES[$step]}"
            job_name="atac_${step_name}_${ds}"
            cpus=$(get_resource ATAC $step CPUS)
            mem=$(get_resource ATAC $step MEM)
            time_limit=$(get_resource ATAC $step TIME)

            echo "  │  Step $step: ${ATAC_STEP_DESC[$step]}"

            submit_step \
                "$job_name" "$ATAC_PARTITION" "$cpus" "$mem" "$time_limit" \
                "$log_dir" "$step_name" "$ATAC_CONDA_ENV_R" \
                "$config_path" "$step" "${REPO_DIR}/atac_only" \
                "${REPO_DIR}/atac_only/slurm/run_step.sbatch" \
                "$prev_job"

            if [ "$RETURNED_JOB_ID" = "FAILED" ]; then
                echo "  │  ✗ Aborting ATAC chain for $ds"
                break
            fi

            atac_chain="${atac_chain:+$atac_chain → }${RETURNED_JOB_ID}"
            prev_job="$RETURNED_JOB_ID"
            TOTAL_JOBS=$((TOTAL_JOBS + 1))
        done

        echo "  │  Chain: ${combined_step2_job:+${combined_step2_job} (combined step 2) → }$atac_chain"
        echo "  └──────────────────────────"
    fi

    echo ""
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
    echo "  Datasets: ${#DATASETS[@]}"
fi

echo ""
echo "  Dependency structure per dataset:"
echo "    combined(1→2→3→4→5→6)"
echo "                ↓ after step 2"
echo "                ├→ rna_only(1→2→3→4)"
echo "                └→ atac_only(1→2→3→4)"
echo ""
echo "  Monitor:  squeue -u $USER"
echo "  Cancel:   scancel <job_id>"
echo "  Logs:     See each pipeline's log directory"
echo ""
