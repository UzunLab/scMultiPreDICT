# Combined (RNA+ATAC) Pipeline - Dataset Configuration Directory

This directory holds **per-dataset configuration files** for multi-dataset runs.

## Quick Start

```bash
# 1. Create a config for each dataset (copy the template)
cp ../R/config_template.R E7.5_rep1.config.R
cp ../R/config_template.R E7.5_rep2.config.R
cp ../R/config_template.R T_Cells.config.R

# 2. Edit each .config.R with dataset-specific paths and parameters
#    Key fields to change per dataset:
#      - SAMPLE_NAME
#      - INPUT_MTX, INPUT_FEATURES, INPUT_BARCODES, INPUT_FRAGMENTS
#      - BASE_OUTPUT_DIR
#      - SPECIES ("mouse" or "human")
#      - DIM_REDUCTION_METHOD
#      - HVG_GENE_FILE / RANDOM_GENE_FILE (or leave empty for auto-selection)

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
- SLURM job names (e.g., `qc_preprocess_E7.5_rep1`)
- Log directory names
- Display in submission output

## Pipeline Steps (serial per dataset)

| Step | Name | Description |
|------|------|-------------|
| 1 | Quality Control | Cell/gene filtering on RNA + ATAC metrics |
| 2 | Data Splitting | Stratified train/val/test + target gene selection |
| 3 | Metacell Creation | kNN smoothing (PCA+LSI, WNN, MultiVI, or scVI+PeakVI) |
| 4 | Feature Extraction | Gene-specific feature matrix (peaks + expression) |
| 5 | Linear Models | Ridge, Lasso, Elastic Net, Random Forest |
| 6 | Neural Network | Feed-forward NN with optional grid search |

## Python Autoencoder Steps (multivi / scvi_peakvi only)

If your `DIM_REDUCTION_METHOD` is `multivi` or `scvi_peakvi`, you need
additional Python preprocessing between steps 2 and 3:

1. Export Seurat to MuData format (03b_export_to_mudata.R)
2. Train autoencoders (python/train_autoencoders.py)

These optional steps are NOT included in the main submit_datasets.sh chain.
Run them separately before starting step 3, or use `--start-step 3` after
completing these manually.

## Example: Different Species Datasets

`E7.5_rep1.config.R` (mouse):
```r
SAMPLE_NAME <- "E7.5_rep1"
SPECIES <- "mouse"
INPUT_MTX <- "/data/E7.5_rep1/matrix.mtx.gz"
INPUT_FRAGMENTS <- "/data/E7.5_rep1/fragments.tsv.gz"
```

`T_Cells.config.R` (human):
```r
SAMPLE_NAME <- "T_Cells"
SPECIES <- "human"
INPUT_MTX <- "/data/T_Cells/matrix.mtx.gz"
INPUT_FRAGMENTS <- "/data/T_Cells/fragments.tsv.gz"
```
