# ============================================================================
# Step_02: Data Splitting for Gene Prediction (Automated)
# ============================================================================
# This script splits preprocessed data into train/validation/test sets.
# It reads configuration from config.R.
#
# Input: Preprocessed Seurat object from Step_01
# Output: Seurat object with split labels + split indices
#
# Usage:
#   1. Ensure config.R is properly configured
#   2. Run: Rscript 02a_data_splitting.R
# ============================================================================

# Source configuration (skip if already loaded by run_pipeline.R)
if (!exists("CONFIG_LOADED")) {
  get_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      return(dirname(normalizePath(sub("^--file=", "", file_arg))))
    }
    return(".")
  }
  script_dir <- get_script_dir()
  
  config_path <- file.path(script_dir, "config.R")
  if (!file.exists(config_path)) {
    config_path <- "config.R"
  }
  if (!file.exists(config_path)) {
    stop("config.R not found! Please copy config_template.R to config.R and edit with your settings.")
  }
  cat("Loading configuration from:", config_path, "\n")
  source(config_path)
}

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# Create output directories
create_output_directories()

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 02: Data Splitting\n")
cat("=", rep("=", 70), "\n\n", sep = "")

# ============================================================================
# LOAD PREPROCESSED DATA
# ============================================================================
cat("=== Loading preprocessed data ===\n\n")

# Input file from Step 01
input_file <- file.path(
  OUTPUT_SEURAT_DIR, 
  paste0(SAMPLE_NAME, "_seurat_multiome_processed.rds")
)

if (!file.exists(input_file)) {
  stop(sprintf("ERROR: Input file not found: %s\nPlease run Step 01 first.", input_file))
}

cat(sprintf("Loading: %s\n", input_file))
seurat_obj <- readRDS(input_file)

cat(sprintf("Loaded object: %d cells, %d genes, %d peaks\n\n", 
            ncol(seurat_obj), 
            nrow(seurat_obj[["RNA"]]),
            nrow(seurat_obj[["ATAC"]])))

# ============================================================================
# DEFINE DATA SPLITS
# ============================================================================
cat("=== Creating train/validation/test splits ===\n\n")

# Get all cell barcodes
all_cells <- colnames(seurat_obj)
n_cells <- length(all_cells)

# Set seed for reproducibility
set.seed(SEED_SPLIT)

# Calculate number of cells for each split
n_train <- floor(SPLIT_TRAIN * n_cells)
n_val <- floor(SPLIT_VAL * n_cells)
n_test <- n_cells - n_train - n_val  # Remainder goes to test

# Shuffle and split
shuffled_idx <- sample(1:n_cells, size = n_cells, replace = FALSE)

train_idx <- shuffled_idx[1:n_train]
val_idx <- shuffled_idx[(n_train + 1):(n_train + n_val)]
test_idx <- shuffled_idx[(n_train + n_val + 1):n_cells]

# Add split labels to Seurat object metadata
seurat_obj$data_split <- NA
seurat_obj$data_split[train_idx] <- "train"
seurat_obj$data_split[val_idx] <- "validation"
seurat_obj$data_split[test_idx] <- "test"

# ============================================================================
# PRINT SPLIT SUMMARY
# ============================================================================
cat("Split configuration:\n")
cat(sprintf("  Target: Train=%.0f%%, Validation=%.0f%%, Test=%.0f%%\n",
            SPLIT_TRAIN * 100, SPLIT_VAL * 100, SPLIT_TEST * 100))
cat(sprintf("  Random seed: %d\n\n", SEED_SPLIT))

cat("Split results:\n")
cat(sprintf("  Train:       %5d cells (%5.1f%%)\n", 
            length(train_idx), 100 * length(train_idx) / n_cells))
cat(sprintf("  Validation:  %5d cells (%5.1f%%)\n", 
            length(val_idx), 100 * length(val_idx) / n_cells))
cat(sprintf("  Test:        %5d cells (%5.1f%%)\n", 
            length(test_idx), 100 * length(test_idx) / n_cells))
cat(sprintf("  Total:       %5d cells\n", n_cells))

# Verify no overlap
all_split_idx <- c(train_idx, val_idx, test_idx)
if (length(unique(all_split_idx)) != n_cells) {
  stop("ERROR: Split indices have overlap or missing cells!")
}
cat("\nVerification: ✓ No overlap between splits\n")

# ============================================================================
# SAVE OUTPUTS
# ============================================================================
cat("\n=== Saving split information ===\n\n")

# 1. Save updated Seurat object with split labels
seurat_output_file <- file.path(
  OUTPUT_SPLITS_DIR, 
  paste0(SAMPLE_NAME, "_seurat_obj_with_splits.rds")
)
saveRDS(seurat_obj, seurat_output_file)
cat(sprintf("Saved: %s\n", seurat_output_file))

# 2. Save split indices as CSV for easy access
split_df <- data.frame(
  cell_barcode = colnames(seurat_obj),
  data_split = seurat_obj$data_split,
  stringsAsFactors = FALSE
)
csv_file <- file.path(OUTPUT_SPLITS_DIR, "cell_data_splits.csv")
write.csv(split_df, csv_file, row.names = FALSE)
cat(sprintf("Saved: %s\n", csv_file))

# 3. Save split indices as RDS for faster loading in R
split_list <- list(
  train = train_idx,
  validation = val_idx,
  test = test_idx,
  config = list(
    sample_name = SAMPLE_NAME,
    split_train = SPLIT_TRAIN,
    split_val = SPLIT_VAL,
    split_test = SPLIT_TEST,
    seed = SEED_SPLIT,
    n_cells = n_cells
  )
)
indices_file <- file.path(OUTPUT_SPLITS_DIR, "split_indices.rds")
saveRDS(split_list, indices_file)
cat(sprintf("Saved: %s\n", indices_file))

# 4. Save split summary as CSV
summary_df <- data.frame(
  split = c("train", "validation", "test", "total"),
  n_cells = c(length(train_idx), length(val_idx), length(test_idx), n_cells),
  percentage = c(
    100 * length(train_idx) / n_cells,
    100 * length(val_idx) / n_cells,
    100 * length(test_idx) / n_cells,
    100
  ),
  stringsAsFactors = FALSE
)
summary_file <- file.path(OUTPUT_SPLITS_DIR, "split_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat(sprintf("Saved: %s\n", summary_file))

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 02 COMPLETE\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Output files:\n")
cat(sprintf("  - %s (Seurat object with split labels)\n", basename(seurat_output_file)))
cat(sprintf("  - %s (cell barcode to split mapping)\n", basename(csv_file)))
cat(sprintf("  - %s (train/val/test indices)\n", basename(indices_file)))
cat(sprintf("  - %s (split statistics)\n", basename(summary_file)))

cat("\nNext step: Run 03a_metacell_creation_pca_lsi.R\n")
cat("           (or train_autoencoder.*.py for non-linear dimensionality reduction)\n")