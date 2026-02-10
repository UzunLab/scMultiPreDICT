# Multi-Dataset Pipeline Guide

This document explains how to run scMultiPreDICT across **multiple datasets in parallel**, with each dataset's steps running **serially** via SLURM dependency chains.

## Architecture Overview

```
For each dataset:   Step 1 → Step 2 → Step 3 → Step 4   (serial, via --dependency=afterok)

Across datasets:    Dataset A: [Step1 → Step2 → Step3 → Step4]
                    Dataset B: [Step1 → Step2 → Step3 → Step4]   (parallel, independent)
                    Dataset C: [Step1 → Step2 → Step3 → Step4]
```

Each dataset gets its own **config file** with dataset-specific paths and parameters. The orchestration script (`submit_datasets.sh`) submits all datasets, using SLURM's `--dependency=afterok:<job_id>` to enforce serial execution within each dataset while allowing all datasets to run independently in parallel.

## Directory Structure (New Files)

```
rna_only/
  datasets/                        # Per-dataset config files
    README.md
    E7.5_rep1.config.R             # ← you create these
    T_Cells.config.R
  slurm/
    run_step.sbatch                # Generic parameterized step runner
  submit_settings.template.sh     # Cluster settings template
  submit_settings.sh              # ← you create (copy from template)
  submit_datasets.sh              # Master orchestration script

atac_only/                         # Same structure
  datasets/
  slurm/
  submit_settings.template.sh
  submit_datasets.sh

combined/                          # Same structure, 6 steps instead of 4
  datasets/
  slurm/
    run_step.sbatch
    run_step_python.sbatch         # For optional Python steps
  submit_settings.template.sh
  submit_datasets.sh
```

## Setup (per pipeline)

### 1. Configure Cluster Settings

```bash
cd rna_only/   # or atac_only/ or combined/
cp submit_settings.template.sh submit_settings.sh
```

Edit `submit_settings.sh`:
- `PARTITION` — your SLURM partition name
- `CONDA_ENV_R` — conda environment with R and Bioconductor
- `LOG_BASE_DIR` — where SLURM logs go
- `STEP_N_CPUS`, `STEP_N_MEM`, `STEP_N_TIME` — per-step resources

### 2. Create Per-Dataset Configs

```bash
# RNA-only / ATAC-only: copy from config_template.R
cp config_template.R datasets/E7.5_rep1.config.R
cp config_template.R datasets/E7.5_rep2.config.R
cp config_template.R datasets/T_Cells.config.R

# Combined: copy from R/config_template.R
cp R/config_template.R datasets/E7.5_rep1.config.R
```

Edit each config with **dataset-specific** values:
- `SAMPLE_NAME` — unique identifier
- Input file paths
- `SPECIES` — `"mouse"` or `"human"`
- Output directories
- Target gene files

### 3. Submit

```bash
chmod +x submit_datasets.sh

# Preview what will happen (no jobs submitted)
./submit_datasets.sh --dry-run

# Submit all datasets
./submit_datasets.sh

# Submit specific datasets only
./submit_datasets.sh --datasets E7.5_rep1,T_Cells

# Start from a specific step (e.g., re-run models only)
./submit_datasets.sh --start-step 3

# Combine options
./submit_datasets.sh --datasets E7.5_rep1 --start-step 2 --stop-step 3
```

## How It Works

1. `submit_datasets.sh` scans `datasets/*.config.R` for dataset configs
2. For each dataset, it submits a chain of SLURM jobs:
   - Step 1 submitted normally
   - Step 2 submitted with `--dependency=afterok:<step1_job_id>`
   - Step 3 submitted with `--dependency=afterok:<step2_job_id>`
   - ...and so on
3. Each SLURM job runs `slurm/run_step.sbatch`, which:
   - Activates the conda environment
   - Runs `Rscript run_pipeline.R --config <dataset_config> --steps <N>`
4. All datasets are independent — if Dataset A fails at step 2, Datasets B and C continue

## Monitoring

```bash
# Check all your running/pending jobs
squeue -u $USER

# Check a specific dataset's jobs
squeue -u $USER -n "metacell_E7.5_rep1"

# View logs
tail -f ~/scMultiPreDICT_logs/rna_only/E7.5_rep1/metacell_12345.log

# Cancel all jobs for a dataset (cancel one, dependents auto-cancel)
scancel <any_job_id_in_chain>

# Cancel ALL your jobs
scancel -u $USER
```

## Pipeline Steps Reference

### Combined Pipeline (6 steps)
| Step | Name | Script |
|------|------|--------|
| 1 | Quality Control | R/01_quality_control.R |
| 2 | Data Splitting | R/02a_data_splitting.R |
| 3 | Metacell Creation | R/03a_metacell_creation_*.R |
| 4 | Feature Extraction | R/04_feature_extraction.R |
| 5 | Linear Models (train + predict) | R/05_linear_tree_models.R |
| 6 | Neural Network (train + predict) | R/06_neural_network.R |

### RNA-only Pipeline (4 steps)
| Step | Name | Script |
|------|------|--------|
| 1 | Metacell Creation | R/03a_metacell_creation.R |
| 2 | Feature Extraction | R/04_feature_extraction.R |
| 3 | Linear Models (train + predict) | R/05_linear_tree_models.R |
| 4 | Neural Network (train + predict) | R/06_neural_network.R |

### ATAC-only Pipeline (4 steps)
| Step | Name | Script |
|------|------|--------|
| 1 | Metacell Creation | R/03a_metacell_creation.R |
| 2 | Feature Extraction | R/04_feature_extraction.R |
| 3 | Linear Models (train + predict) | R/05_linear_tree_models.R |
| 4 | Neural Network (train + predict) | R/06_neural_network.R |

## Backwards Compatibility

The original single-dataset workflow still works:
- Edit `config.R` directly
- Run `Rscript run_pipeline.R`

The new multi-dataset system is additive — it doesn't change any existing files.

## Cross-Pipeline Orchestration (submit_all.sh)

The RNA-only and ATAC-only pipelines depend on the **combined pipeline's step 2** output (Seurat splits and target gene lists). The `submit_all.sh` script at the repo root handles this dependency automatically:

```
COMBINED:  step1 → step2 → step3 → step4 → step5 → step6
                      ↓ (after step 2 completes)
RNA_ONLY:           rna1 → rna2 → rna3 → rna4
                      ↓ (after step 2 completes)
ATAC_ONLY:         atac1 → atac2 → atac3 → atac4
```

After combined step 2 finishes, three branches run **in parallel**: combined steps 3-6, RNA-only steps 1-4, and ATAC-only steps 1-4.

### Setup

1. Create `submit_settings.sh` for **each pipeline** you want to run
2. Create matching dataset configs in each pipeline's `datasets/` directory
3. Run from the repo root:

```bash
./submit_all.sh --dry-run                  # Preview the full graph
./submit_all.sh                            # Submit everything
./submit_all.sh --pipelines combined,rna   # Skip ATAC
./submit_all.sh --datasets T_Cells         # Specific datasets only
```

**Important:** If a dataset has no config in `rna_only/datasets/` or `atac_only/datasets/`, that pipeline is automatically skipped for that dataset.

## Troubleshooting

**"No .config.R files found"**: Create at least one config in the `datasets/` directory.

**"Settings file not found"**: Copy `submit_settings.template.sh` to `submit_settings.sh`.

**Job fails, rest of chain stays pending**: SLURM's `afterok` dependency means downstream jobs won't start if a predecessor fails. Fix the issue and resubmit from the failed step:
```bash
./submit_datasets.sh --datasets FAILED_DATASET --start-step FAILED_STEP_NUM
```

**Different resources per dataset**: Resource allocation is per-step (not per-dataset). If a specific dataset needs more memory, adjust `submit_settings.sh` to accommodate the largest dataset, or submit that dataset separately with modified settings.
