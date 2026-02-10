# ATAC-only Pipeline - Dataset Configuration Directory

This directory holds **per-dataset configuration files** for multi-dataset runs.

## Quick Start

```bash
# 1. Create a config for each dataset (copy the template)
cp ../config_template.R E7.5_rep1.config.R
cp ../config_template.R E7.5_rep2.config.R
cp ../config_template.R T_Cells.config.R

# 2. Edit each .config.R with dataset-specific paths and parameters
#    Key fields to change per dataset:
#      - SAMPLE_NAME
#      - INPUT_SEURAT_SPLITS (path to Seurat .rds from combined pipeline)
#      - INPUT_TARGET_GENES_DIR
#      - BASE_OUTPUT_DIR
#      - SPECIES ("mouse" or "human")
#      - HVG_GENE_FILE / RANDOM_GENE_FILE

# 3. Configure cluster settings
cd ..
cp submit_settings.template.sh submit_settings.sh
# Edit submit_settings.sh with your partition, conda env, etc.

# 4. Submit all datasets
./submit_datasets.sh              # all datasets
./submit_datasets.sh --dry-run    # preview first
```

## File Naming Convention

Each config file must be named: `<DATASET_NAME>.config.R`

The dataset name (filename without `.config.R`) is used for:
- SLURM job names (e.g., `metacell_E7.5_rep1`)
- Log directory names
- Display in submission output

## Important Notes

- The ATAC-only pipeline requires outputs from the **combined pipeline**
  (specifically: Seurat object with splits and target gene lists).
  Run combined pipeline steps 1-2 for each dataset first.
- Each dataset runs its 4 steps serially (with SLURM dependencies).
- Different datasets run in parallel (no cross-dataset dependencies).
