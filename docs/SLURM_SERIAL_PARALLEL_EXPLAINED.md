# Understanding Serial & Parallel SLURM Pipelines — A Beginner's Guide

This document explains how the multi-dataset submission system works in scMultiPreDICT, starting from the SLURM basics.

---

## Table of Contents

1. [The Problem](#the-problem)
2. [Two Key SLURM Concepts](#two-key-slurm-concepts)
   - [Serial Execution (Dependency Chains)](#1-serial-execution-dependency-chains)
   - [Parallel Execution (Independent Chains)](#2-parallel-execution-independent-chains)
3. [What Was Built](#what-was-built)
   - [Piece 1: Generic Sbatch Template](#piece-1-one-generic-sbatch-template-slurmrun_stepsbatch)
   - [Piece 2: Per-Dataset Config Files](#piece-2-per-dataset-config-files-datasetsconfigr)
   - [Piece 3: The Orchestrator Script](#piece-3-the-orchestrator-submit_datasetssh)
   - [Cluster Settings](#cluster-settings-submit_settingssh)
4. [The Full Flow](#the-full-flow-when-you-run-it)
5. [How to Add a New Dataset](#how-to-add-a-new-dataset-in-the-future)
6. [Useful Flags](#useful-flags-to-remember)
7. [Monitoring & Recovery](#monitoring--recovery)
8. [File Reference](#file-reference)

---

## The Problem

The pipeline has **steps that must run in order** (step 1 finishes before step 2 starts, etc.) because each step's output is the next step's input. Before, you'd manually submit each SLURM job one at a time:

```bash
sbatch step_030.sbatch    # wait for it to finish...
sbatch step_040.sbatch    # wait again...
sbatch step_050.sbatch    # wait again...
sbatch step_060.sbatch
```

And this only handled **one dataset**. If you had 3 datasets, you'd have to do this 3 times, waiting around.

---

## Two Key SLURM Concepts

### 1. Serial Execution (Dependency Chains)

When you submit a job with `sbatch`, SLURM gives you back a **job ID** (e.g., `12345`). You can tell the *next* job: "don't start until job 12345 succeeds" using the `--dependency` flag:

```bash
# Submit step 1 — runs immediately
JOB1=$(sbatch --parsable step1.sbatch)
# JOB1 now contains something like "12345"

# Submit step 2 — waits for step 1 to finish successfully
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 step2.sbatch)
# "afterok:12345" means "only start after job 12345 exits with code 0 (success)"

# Submit step 3 — waits for step 2
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 step3.sbatch)
```

This creates a **chain**: `Step1 → Step2 → Step3`. SLURM manages the waiting — you submit all three at once and walk away. The jobs sit in "pending" state with reason `(Dependency)` until their predecessor finishes.

**If step 1 fails**, steps 2 and 3 will **never run** (they stay pending forever until you cancel them). This is exactly what you want — no point extracting features if metacell creation crashed.

#### Key flag reference

| Flag | Meaning |
|------|---------|
| `--dependency=afterok:JOB_ID` | Start only after JOB_ID completes **successfully** (exit code 0) |
| `--dependency=afterany:JOB_ID` | Start after JOB_ID finishes, **regardless** of success/failure |
| `--dependency=after:JOB_ID` | Start after JOB_ID **begins** (not finishes) |

We use `afterok` because each step needs the previous step to have **succeeded**.

### 2. Parallel Execution (Independent Chains)

If you have two datasets, their pipelines have **no relationship** to each other. Dataset A's metacell step doesn't need to wait for Dataset B's metacell step. So you just submit two independent chains:

```bash
# Dataset A chain
A1=$(sbatch --parsable step1_A.sbatch)
A2=$(sbatch --parsable --dependency=afterok:$A1 step2_A.sbatch)
A3=$(sbatch --parsable --dependency=afterok:$A2 step3_A.sbatch)

# Dataset B chain (no dependency on anything from A!)
B1=$(sbatch --parsable step1_B.sbatch)
B2=$(sbatch --parsable --dependency=afterok:$B1 step2_B.sbatch)
B3=$(sbatch --parsable --dependency=afterok:$B2 step3_B.sbatch)
```

SLURM schedules everything optimally. If the cluster has room, both A1 and B1 run at the same time. The result looks like:

```
Time →

Dataset A:  [==Step1==]  [===Step2===]  [==Step3==]
Dataset B:  [==Step1==]  [===Step2===]  [==Step3==]
                ↑ both run simultaneously (parallel across datasets)
                  ↑ each dataset's steps wait for the previous (serial within)
```

**That's it.** Serial = `--dependency=afterok:PREV_ID`. Parallel = no dependency between them.

---

## What Was Built

### The problem with the old approach

The old sbatch files (like `step_030_YOUR_SAMPLE_NAME.sbatch`) had everything **hardcoded** — the sample name, paths, conda env, resource amounts. For each new dataset, you'd need to copy and edit 4 sbatch files. For 5 datasets × 4 steps = 20 files to manage. That doesn't scale.

### The solution: 3 pieces

```
┌─────────────────────────────────────────────────────────────────────┐
│                       submit_datasets.sh                            │
│                      (the orchestrator)                             │
│                                                                     │
│  Reads submit_settings.sh  →  cluster config (partition, RAM, etc.) │
│  Reads datasets/*.config.R →  one per dataset                      │
│  Submits slurm/run_step.sbatch with --dependency chains            │
└─────────────────────────────────────────────────────────────────────┘
```

#### Piece 1: One generic sbatch template (`slurm/run_step.sbatch`)

Instead of 4 hardcoded sbatch files per dataset, there's now **one template** that accepts **environment variables**:

```bash
# Old way — hardcoded:
#SBATCH --job-name=metacell_E7.5_rep1
Rscript 03a_metacell_creation.R

# New way — parameterized:
# The job receives CONFIG_PATH and STEP_NUM as environment variables
Rscript run_pipeline.R --config "$CONFIG_PATH" --steps "$STEP_NUM"
```

The key insight: your `run_pipeline.R` already supports `--config PATH --steps N`. So the sbatch template just calls that with the variables passed in. You never need to edit this file.

**How variables get passed in:** The `--export` flag on `sbatch` sets environment variables inside the job:

```bash
sbatch --export="ALL,CONFIG_PATH=/path/to/E7.5_rep1.config.R,STEP_NUM=1,PIPELINE_DIR=/path/to/rna_only,CONDA_ENV_NAME=r-bioc-43" slurm/run_step.sbatch
```

Inside `run_step.sbatch`, `$CONFIG_PATH`, `$STEP_NUM`, etc. are available as regular bash variables.

#### Piece 2: Per-dataset config files (`datasets/*.config.R`)

Instead of one `config.R` that you keep overwriting, each dataset gets its **own** config file:

```
datasets/
  E7.5_rep1.config.R    ← paths/settings for E7.5_rep1
  E7.5_rep2.config.R    ← paths/settings for E7.5_rep2
  T_Cells.config.R      ← paths/settings for T_Cells
```

Each is a copy of `config_template.R` with dataset-specific values filled in (SAMPLE_NAME, file paths, species, etc.). They live side by side and never conflict.

#### Piece 3: The orchestrator (`submit_datasets.sh`)

This bash script automates the "submit chains" pattern. Here's the logic in pseudocode:

```
load cluster settings from submit_settings.sh
find all *.config.R files in datasets/

for each config file:
    dataset_name = filename without .config.R
    prev_job_id = ""

    for each step (1, 2, 3, 4):
        if prev_job_id exists:
            job_id = sbatch --dependency=afterok:prev_job_id ...
        else:
            job_id = sbatch ...  (first step, no dependency)

        prev_job_id = job_id

    print "Chain: job1 → job2 → job3 → job4"
```

Because the outer loop processes each dataset independently, dataset chains have no cross-dependencies and run **in parallel**.

#### Cluster settings (`submit_settings.sh`)

This separates **cluster-specific stuff** (partition name, conda env, how much RAM per step) from **dataset-specific stuff** (paths, species, etc.):

```bash
# submit_settings.sh — configured ONCE for your cluster
PARTITION="compute"
CONDA_ENV_R="r-bioc-43"
STEP_1_CPUS=32
STEP_1_MEM="128G"
STEP_1_TIME="06:00:00"
# ...
```

You configure cluster settings once, dataset configs once each, and you're done.

---

## The Full Flow When You Run It

```bash
./submit_datasets.sh
```

1. Reads `submit_settings.sh` for cluster config (partition, RAM, etc.)
2. Scans `datasets/*.config.R` — finds 3 files → 3 datasets
3. For **E7.5_rep1**: submits step 1 (gets job 100), step 2 depends on 100 (gets 101), step 3 depends on 101 (gets 102), step 4 depends on 102 (gets 103)
4. For **E7.5_rep2**: submits step 1 (gets 104), step 2 depends on 104 (gets 105), ... (**no link to jobs 100–103!**)
5. For **T_Cells**: same pattern, independent chain
6. Prints a summary showing all the chains

From SLURM's perspective, 12 jobs are now queued:

```
Job 100: metacell_E7.5_rep1         (running)
Job 101: feature_extract_E7.5_rep1  (pending — Dependency on 100)
Job 102: linear_models_E7.5_rep1    (pending — Dependency on 101)
Job 103: neural_network_E7.5_rep1   (pending — Dependency on 102)
Job 104: metacell_E7.5_rep2         (running)       ← parallel with 100!
Job 105: feature_extract_E7.5_rep2  (pending — Dependency on 104)
...
```

SLURM runs as many as the cluster allows simultaneously. Each dataset's steps wait for their predecessor, but different datasets run in parallel.

---

## How to Add a New Dataset in the Future

```bash
cd rna_only/   # or atac_only/ or combined/

# 1. Copy the template
cp config_template.R datasets/NEW_SAMPLE.config.R

# 2. Edit it — change SAMPLE_NAME, file paths, species, etc.
nano datasets/NEW_SAMPLE.config.R    # or use any editor

# 3. Submit (it automatically picks up the new file)
./submit_datasets.sh --dry-run    # preview first!
./submit_datasets.sh              # for real
```

That's it. No sbatch files to create, no dependency chains to manage manually.

---

## Useful Flags to Remember

```bash
# Preview what would happen (NO jobs submitted)
./submit_datasets.sh --dry-run

# Submit ALL datasets in datasets/
./submit_datasets.sh

# Submit only specific datasets
./submit_datasets.sh --datasets E7.5_rep1
./submit_datasets.sh --datasets E7.5_rep1,T_Cells

# Start from a specific step (e.g., steps 1-2 already done, rerun from 3)
./submit_datasets.sh --start-step 3

# Stop after a specific step
./submit_datasets.sh --stop-step 2

# Run ONLY step 3
./submit_datasets.sh --start-step 3 --stop-step 3

# Combine everything
./submit_datasets.sh --datasets E7.5_rep1 --start-step 2 --stop-step 3 --dry-run
```

---

## Monitoring & Recovery

### Check job status

```bash
# See all your running/pending jobs
squeue -u $USER

# Example output:
#   JOBID   NAME                       STATE     REASON
#   100     metacell_E7.5_rep1         RUNNING   None
#   101     feature_extract_E7.5_rep1  PENDING   (Dependency)
#   104     metacell_E7.5_rep2         RUNNING   None
#   105     feature_extract_E7.5_rep2  PENDING   (Dependency)
```

### View logs

```bash
# Logs are at: LOG_BASE_DIR/<dataset_name>/<step_name>_<job_id>.log
tail -f ~/scMultiPreDICT_logs/rna_only/E7.5_rep1/metacell_100.log

# Error output:
cat ~/scMultiPreDICT_logs/rna_only/E7.5_rep1/metacell_100.err
```

### If a step fails

When a step fails (non-zero exit code), `afterok` dependencies mean **all downstream jobs for that dataset stay pending forever**. Other datasets are unaffected.

```bash
# 1. Check what went wrong
cat ~/scMultiPreDICT_logs/rna_only/E7.5_rep1/metacell_100.err

# 2. Cancel the stuck pending jobs in that chain
scancel 101 102 103
# Or cancel ALL your jobs: scancel -u $USER

# 3. Fix the issue (edit config, fix data, etc.)

# 4. Resubmit from the failed step
./submit_datasets.sh --datasets E7.5_rep1 --start-step 1
```

### Cancel jobs

```bash
scancel 12345              # cancel a specific job
scancel -u $USER           # cancel ALL your jobs
scancel -n metacell_E7.5_rep1  # cancel by job name
```

When you cancel a job, its dependent jobs also get cancelled automatically.

---

## File Reference

### New files added per pipeline (`rna_only/`, `atac_only/`, `combined/`)

```
<pipeline>/
│
├── submit_datasets.sh              # Master orchestration script (run this!)
├── submit_settings.template.sh     # Cluster settings template
├── submit_settings.sh              # YOUR cluster settings (you create this)
│
├── datasets/                       # Per-dataset configuration directory
│   ├── README.md                   # Setup instructions
│   ├── .gitkeep                    # Keeps empty dir in git
│   ├── E7.5_rep1.config.R          # (you create these)
│   ├── E7.5_rep2.config.R
│   └── T_Cells.config.R
│
├── slurm/
│   ├── run_step.sbatch             # Generic parameterized step runner
│   └── run_step_python.sbatch      # (combined/ only) For Python steps
│
├── config_template.R               # Template to copy into datasets/
├── config.R                        # Still works for single-dataset runs
├── run_pipeline.R                  # Unchanged — core pipeline runner
└── R/                              # Unchanged — analysis scripts
```

### What each file does

| File | You edit it? | Purpose |
|------|:---:|---------|
| `submit_settings.template.sh` | No | Reference template |
| `submit_settings.sh` | **Yes, once** | Your cluster's partition, conda env, RAM per step |
| `datasets/SAMPLE.config.R` | **Yes, per dataset** | Sample name, input paths, species, output dirs |
| `submit_datasets.sh` | No | Reads settings + configs, submits SLURM chains |
| `slurm/run_step.sbatch` | No | Receives env vars, runs `Rscript run_pipeline.R --config X --steps N` |

### Pipeline step counts

| Pipeline | Steps | What they are |
|----------|:-----:|---------------|
| Combined (RNA+ATAC) | 6 | QC → Split → Metacell → Features → Linear Models → Neural Net |
| RNA-only | 4 | Metacell → Features → Linear Models → Neural Net |
| ATAC-only | 4 | Metacell → Features → Linear Models → Neural Net |

---

## Summary

The core pattern is simple:

- **Serial within a dataset:** `--dependency=afterok:PREVIOUS_JOB_ID`
- **Parallel across datasets:** Submit independent chains with no cross-dependencies
- **Parameterized jobs:** One sbatch template + `--export` variables instead of many hardcoded files
- **Per-dataset configs:** Each dataset gets its own `.config.R` in the `datasets/` directory

This pattern applies to **any** multi-step SLURM workflow — not just scMultiPreDICT. Once you understand `--dependency` and `--export`, you can build this for any pipeline.
