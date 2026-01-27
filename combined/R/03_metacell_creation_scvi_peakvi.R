# ============================================================================
# Step_030: Metacell Creation via KNN Smoothing - scVI + PeakVI Latent (Automated)
# ============================================================================
# This script creates smoothed expression profiles using k-NN smoothing
# in a joint latent space from scVI (RNA) + PeakVI (ATAC). 
# It reads configuration from config.R.
#
# APPROACH:
# - Load pre-computed scVI (RNA) and PeakVI (ATAC) latent embeddings
# - Z-score standardize each latent space separately
# - Concatenate to form joint space
# - For each cell, find k nearest neighbors in joint space
# - Smooth expression by averaging with k neighbors (preserves all cells)
# - Training: fit z-score parameters
# - Val/Test: apply training z-score, smooth among val/test cells only
#
# PREREQUISITES:
# - Run train_autoeencoder.py to generate scVI and PeakVI latent embeddings
# - Latent files: latent_scvi_rna_*.csv, latent_peakvi_atac_*.csv
#
# Input: Seurat object with split labels from Step_020
# Output: Smoothed expression matrices for train/val/test splits
#
# Usage:
#   1. Ensure config.R is properly configured with SCVI_LATENT_DIR/PEAKVI_LATENT_DIR
#   2. Run: Rscript 03_metacell_creation_scvi_peakvi.R
# ============================================================================

# Get the directory of this script to source config.R
# Use commandArgs() which works reliably with Rscript
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  # Fallback for interactive use or source()
  return(".")
}
script_dir <- get_script_dir()

# Source configuration file
config_path <- file.path(script_dir, "config_template.R")
if (!file.exists(config_path)) {
  config_path <- "config_template.R"
}

if (!file.exists(config_path)) {
  stop("config_template.R not found! Please ensure config_template.R is in the same directory as this script.")
}

cat("Loading configuration from:", config_path, "\n")
source(config_path)

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(RANN)
  library(dplyr)
  library(cowplot)
  library(reshape2)
})

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

#' Plotting ggplot2 theme
theme_pub <- function(base_size = 12, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5, 
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "grey30",
                                   margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 2, hjust = 1, color = "grey50"),
      axis.title = element_text(face = "bold", size = base_size),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(face = "bold", size = base_size, color = "black"),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      plot.margin = margin(15, 15, 15, 15),
      plot.tag = element_text(face = "bold", size = base_size + 4)
    )
}

#' Colorblind-friendly palette
colorblind_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

#' Split-specific colors (consistent across all plots)
split_colors_pub <- c(
  "train" = "#0072B2",       # Blue
  "validation" = "#E69F00",  # Orange  
  "test" = "#009E73"         # Green
)

#' Modality colors
modality_colors_pub <- c(
  "scVI (RNA)" = "#56B4E9", 
  "PeakVI (ATAC)" = "#2C786C"
)

#' Save plot in both PNG and PDF formats
save__plot <- function(plot, filename, width = 7, height = 5, dpi = 600) {
  ggsave(paste0(filename, ".png"), plot = plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(filename, ".pdf"), plot = plot, width = width, height = height, device = cairo_pdf)
  invisible(plot)
}

# Create output directories
create_output_directories()

# Set parameters from config
k_neighbors <- K_NEIGHBORS
seed <- SEED_METACELL

set.seed(seed)

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 050: Metacell Creation (scVI + PeakVI k-NN Smoothing)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Parameters:\n")
cat(sprintf("  k-NN neighbors: %d\n", k_neighbors))
cat(sprintf("  Random seed: %d\n\n", seed))

# ============================================================================
# DEFINE LATENT SPACE DIRECTORY
# ============================================================================
# scVI and PeakVI latent files should be in this directory
if (exists("SCVI_PEAKVI_LATENT_DIR") && SCVI_PEAKVI_LATENT_DIR != "") {
  latent_dir <- path.expand(SCVI_PEAKVI_LATENT_DIR)
} else {
  # Default location based on output structure
  latent_dir <- file.path(OUTPUT_LATENT_DIR)
}

cat("scVI/PeakVI latent directory:", latent_dir, "\n\n")

# ============================================================================
# LOAD PREPROCESSED DATA
# ============================================================================
cat("=== Loading preprocessed data ===\n\n")

input_file <- file.path(
  OUTPUT_SPLITS_DIR, 
  paste0(SAMPLE_NAME, "_seurat_obj_with_splits.rds")
)

if (!file.exists(input_file)) {
  stop(sprintf("ERROR: Input file not found: %s\nPlease run Step_020 first.", input_file))
}

cat(sprintf("Loading: %s\n", input_file))
seurat_obj <- readRDS(input_file)

cat(sprintf("Loaded object: %d cells, %d genes, %d peaks\n\n", 
            ncol(seurat_obj), 
            nrow(seurat_obj[["RNA"]]),
            nrow(seurat_obj[["ATAC"]])))

# ============================================================================
# VALIDATE INPUTS
# ============================================================================
cat("=== Validating inputs ===\n\n")

if (!"RNA" %in% names(seurat_obj@assays)) {
  stop("ERROR: RNA assay not found in Seurat object")
}
if (!"ATAC" %in% names(seurat_obj@assays)) {
  stop("ERROR: ATAC assay not found in Seurat object")
}
if (!"data_split" %in% colnames(seurat_obj@meta.data)) {
  stop("ERROR: data_split column not found. Please run Step_020 first.")
}

# Check cells per split
for (split in c("train", "validation", "test")) {
  n_cells <- sum(seurat_obj$data_split == split)
  cat(sprintf("  %s: %d cells\n", split, n_cells))
  if (n_cells == 0) {
    stop(sprintf("ERROR: No cells found for %s split", split))
  }
  if (n_cells < k_neighbors) {
    warning(sprintf("WARNING: %s has only %d cells (< k=%d). k will be adjusted.", 
                    split, n_cells, k_neighbors))
  }
}
cat("\nValidation passed ✓\n")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Z-score functions
z_fit_cols <- function(M) { 
  mu <- colMeans(M)
  sd <- apply(M, 2, sd)
  sd[!is.finite(sd) | sd < 1e-8] <- 1
  list(mu = mu, sd = sd)
}

z_apply_cols <- function(M, mu, sd) { 
  sweep(sweep(M, 2, mu, "-"), 2, sd, "/")
}

# Build k-NN smoothing matrix
build_knn_smooth <- function(scores, cell_names, k_smooth = 20, seed = 1) {
  set.seed(seed)
  n <- nrow(scores)
  total_k <- k_smooth
  k_use <- max(1, min(total_k - 1, n - 1))
  
  if (n <= 1) {
    warning(sprintf("Too few cells (%d) for k-NN smoothing. Using identity.", n))
    return(sparseMatrix(i = 1, j = 1, x = 1, dims = c(n, n), dimnames = list(cell_names, cell_names)))
  }
  
  nn <- RANN::nn2(data = scores, query = scores, k = k_use + 1, searchtype = "standard")
  idx <- nn$nn.idx
  total_k_eff <- k_use + 1
  weight  <- 1 / total_k_eff
  
  i_idx <- as.vector(t(idx))
  j_idx <- rep(seq_len(n), each = total_k)
  x_vals <- rep(weight, length(i_idx))
  
  A <- sparseMatrix(i = i_idx, j = j_idx, x = x_vals, dims = c(n, n), dimnames = list(cell_names, cell_names))
  
  stopifnot(identical(rownames(A), cell_names), identical(colnames(A), cell_names))
  A
}

# Apply k-NN smoothing to expression
smooth_expression <- function(rna_counts, atac_counts, A_smooth) {
  rna_smoothed <- rna_counts %*% A_smooth
  atac_smoothed <- atac_counts %*% A_smooth
  
  rna_lib <- Matrix::colSums(rna_smoothed)
  rna_cpm <- t(t(rna_smoothed) / (rna_lib + 1e-8)) * 1e6
  rna_log1p <- log1p(rna_cpm)
  
  atac_lib <- Matrix::colSums(atac_smoothed)
  atac_cpm <- t(t(atac_smoothed) / (atac_lib + 1e-8)) * 1e6
  atac_log1p <- log1p(atac_cpm)
  
  list(
    rna_log1p = rna_log1p,
    rna_lib = rna_lib,
    atac_counts = atac_smoothed,
    atac_lib = atac_lib,
    atac_log1p = atac_log1p
  )
}

# ============================================================================
# PROCESS TRAINING DATA
# ============================================================================
cat("\n=== Processing TRAINING data ===\n\n")

train_cells <- colnames(seurat_obj)[seurat_obj$data_split == "train"]
seurat_train <- seurat_obj[, train_cells]

cat(sprintf("Training cells: %d\n", length(train_cells)))

# Ensure RNA counts are sparse
DefaultAssay(seurat_train) <- "RNA"
rc <- GetAssayData(seurat_train, assay = "RNA", layer = "counts")
if (!inherits(rc, "dgCMatrix")) {
  seurat_train <- SetAssayData(seurat_train, assay = "RNA", layer = "counts", new.data = as(rc, "dgCMatrix"))
}

cat(sprintf("Genes: %d, Peaks: %d\n", nrow(seurat_train[["RNA"]]), nrow(seurat_train[["ATAC"]])))

# Load scVI latent for TRAINING (RNA)
cat("\n[INFO] Loading scVI (RNA) latent for training...\n")
scvi_train_file <- file.path(latent_dir, "latent_scvi_rna_train.csv")

if (!file.exists(scvi_train_file)) {
  stop(sprintf("ERROR: scVI latent file not found: %s\nPlease run Step_045 (Python scVI) first.", scvi_train_file))
}

latent_rna_train <- read.csv(scvi_train_file, row.names = 1, check.names = FALSE)

# Keep only the scVI latent columns
rna_scores_train <- as.matrix(latent_rna_train[, grepl("^scVI_RNA_", colnames(latent_rna_train)), drop = FALSE])
rna_scores_train <- rna_scores_train[train_cells, , drop = FALSE]

cat(sprintf("scVI RNA latent (train): %d cells × %d dims\n",
            nrow(rna_scores_train), ncol(rna_scores_train)))

# Load PeakVI latent for TRAINING (ATAC)
cat("[INFO] Loading PeakVI (ATAC) latent for training...\n")
peakvi_train_file <- file.path(latent_dir, "latent_peakvi_atac_train.csv")

if (!file.exists(peakvi_train_file)) {
  stop(sprintf("ERROR: PeakVI latent file not found: %s\nPlease run Step_045 (Python PeakVI) first.", peakvi_train_file))
}

latent_atac_train <- read.csv(peakvi_train_file, row.names = 1, check.names = FALSE)

atac_scores_train <- as.matrix(latent_atac_train[, grepl("^PeakVI_ATAC_", colnames(latent_atac_train)), drop = FALSE])
atac_scores_train <- atac_scores_train[train_cells, , drop = FALSE]

cat(sprintf("PeakVI ATAC latent (train): %d cells × %d dims\n",
            nrow(atac_scores_train), ncol(atac_scores_train)))

# Z-score standardization (fit on training data)
cat("\n[INFO] Computing z-score parameters from training...\n")
rna_zs <- z_fit_cols(rna_scores_train)
atac_zs <- z_fit_cols(atac_scores_train)

rna_train_z <- z_apply_cols(rna_scores_train, rna_zs$mu, rna_zs$sd)
atac_train_z <- z_apply_cols(atac_scores_train, atac_zs$mu, atac_zs$sd)

# Create joint space by concatenating
joint_train <- cbind(rna_train_z, atac_train_z)

cat(sprintf("Joint space (scVI + PeakVI): %d cells × %d dimensions\n", 
            nrow(joint_train), ncol(joint_train)))

# Build k-NN smoothing matrix for TRAINING
cat(sprintf("\n[INFO] Building k-NN smoothing matrix (k=%d) for training...\n", k_neighbors))
A_smooth_train <- build_knn_smooth(joint_train, train_cells, k_smooth = k_neighbors, seed = seed)

neighbors_per_cell <- Matrix::rowSums(A_smooth_train > 0)
cat(sprintf("k-NN stats: mean=%.1f, min=%d, max=%d neighbors per cell\n",
            mean(neighbors_per_cell), min(neighbors_per_cell), max(neighbors_per_cell)))

# Smooth TRAINING expression
cat("\n[INFO] Smoothing training expression...\n")
rna_counts_train <- GetAssayData(seurat_train, assay = "RNA", layer = "counts")
atac_counts_train <- GetAssayData(seurat_train, assay = "ATAC", layer = "counts")

smoothed_train <- smooth_expression(rna_counts_train, atac_counts_train, A_smooth_train)
cat(sprintf("[OK] Training: %d cells smoothed (preserved all cells)\n", ncol(smoothed_train$rna_log1p)))

# Save training transformation parameters
train_transforms <- list(
  rna_zs = rna_zs,
  atac_zs = atac_zs,
  k_neighbors = k_neighbors,
  rna_latent_dims = ncol(rna_scores_train),
  atac_latent_dims = ncol(atac_scores_train)
)
saveRDS(train_transforms, file.path(OUTPUT_METACELLS_DIR, "train_transforms_scvi_peakvi.rds"))
cat("[OK] Saved train_transforms_scvi_peakvi.rds\n")

# Save smoothed training data
train_output <- list(
  rna_log1p = smoothed_train$rna_log1p,
  rna_lib = smoothed_train$rna_lib,
  atac_counts = smoothed_train$atac_counts,
  atac_lib = smoothed_train$atac_lib,
  atac_log1p = smoothed_train$atac_log1p,
  rna_scvi_latent = rna_scores_train,
  atac_peakvi_latent = atac_scores_train,
  joint_space = joint_train,
  cells = train_cells,
  split_name = "train",
  processing_seed = seed,
  k_smooth = k_neighbors,
  method = "scVI_PeakVI"
)
saveRDS(train_output, file.path(OUTPUT_METACELLS_DIR, "smoothed_train.rds"))
cat("[OK] Saved smoothed_train.rds\n")

# ============================================================================
# PROCESS VALIDATION AND TEST DATA
# ============================================================================
cat("\n=== Processing VALIDATION and TEST data ===\n")

train_transforms <- readRDS(file.path(OUTPUT_METACELLS_DIR, "train_transforms_scvi_peakvi.rds"))

for (split_name in c("validation", "test")) {
  cat(sprintf("\n--- Processing %s ---\n", toupper(split_name)))
  
  split_cells <- colnames(seurat_obj)[seurat_obj$data_split == split_name]
  seurat_split <- seurat_obj[, split_cells]
  
  cat(sprintf("%s cells: %d\n", split_name, length(split_cells)))
  
  # Ensure data is sparse
  DefaultAssay(seurat_split) <- "RNA"
  
  rc <- GetAssayData(seurat_split, assay = "RNA", layer = "counts")
  if (!inherits(rc, "dgCMatrix")) {
    seurat_split <- SetAssayData(seurat_split, assay = "RNA", layer = "counts", new.data = as(rc, "dgCMatrix"))
  }
  
  ac <- GetAssayData(seurat_split, assay = "ATAC", layer = "counts")
  if (!inherits(ac, "dgCMatrix")) {
    seurat_split <- SetAssayData(seurat_split, assay = "ATAC", layer = "counts", new.data = as(ac, "dgCMatrix"))
  }
  
  # Load scVI / PeakVI latent for this split
  if (split_name == "validation") {
    rna_file <- file.path(latent_dir, "latent_scvi_rna_val.csv")
    atac_file <- file.path(latent_dir, "latent_peakvi_atac_val.csv")
  } else {
    rna_file <- file.path(latent_dir, "latent_scvi_rna_test.csv")
    atac_file <- file.path(latent_dir, "latent_peakvi_atac_test.csv")
  }
  
  cat(sprintf("[INFO] Loading %s scVI RNA latent...\n", split_name))
  if (!file.exists(rna_file)) {
    stop(sprintf("ERROR: scVI latent file not found: %s", rna_file))
  }
  latent_rna_split <- read.csv(rna_file, row.names = 1, check.names = FALSE)
  rna_scores_split <- as.matrix(latent_rna_split[, grepl("^scVI_RNA_", colnames(latent_rna_split)), drop = FALSE])
  rna_scores_split <- rna_scores_split[split_cells, , drop = FALSE]
  
  cat(sprintf("[INFO] Loading %s PeakVI ATAC latent...\n", split_name))
  if (!file.exists(atac_file)) {
    stop(sprintf("ERROR: PeakVI latent file not found: %s", atac_file))
  }
  latent_atac_split <- read.csv(atac_file, row.names = 1, check.names = FALSE)
  atac_scores_split <- as.matrix(latent_atac_split[, grepl("^PeakVI_ATAC_", colnames(latent_atac_split)), drop = FALSE])
  atac_scores_split <- atac_scores_split[split_cells, , drop = FALSE]
  
  cat(sprintf("scVI RNA latent (%s): %d cells × %d dims\n",
              split_name, nrow(rna_scores_split), ncol(rna_scores_split)))
  cat(sprintf("PeakVI ATAC latent (%s): %d cells × %d dims\n",
              split_name, nrow(atac_scores_split), ncol(atac_scores_split)))
  
  # Apply TRAINING z-score parameters
  cat("[INFO] Applying training z-score normalization...\n")
  rna_split_z <- z_apply_cols(rna_scores_split, train_transforms$rna_zs$mu, train_transforms$rna_zs$sd)
  atac_split_z <- z_apply_cols(atac_scores_split, train_transforms$atac_zs$mu, train_transforms$atac_zs$sd)
  
  # Create joint space (in training coordinate system)
  joint_split <- cbind(rna_split_z, atac_split_z)
  
  cat(sprintf("Joint space (%s): %d cells × %d dimensions\n",
              split_name, nrow(joint_split), ncol(joint_split)))
  
  # Build k-NN smoothing matrix WITHIN this split
  cat(sprintf("[INFO] Building k-NN within %s cells (k=%d, no training data)...\n", split_name, k_neighbors))
  A_smooth_split <- build_knn_smooth(joint_split, split_cells, k_smooth = k_neighbors, seed = seed)
  
  neighbors_per_cell <- Matrix::rowSums(A_smooth_split > 0)
  cat(sprintf("k-NN stats: mean=%.1f, min=%d, max=%d neighbors\n",
              mean(neighbors_per_cell), min(neighbors_per_cell), max(neighbors_per_cell)))
  
  # Smooth using ONLY cells from this split
  cat(sprintf("[INFO] Smoothing %s expression...\n", split_name))
  rna_counts_split <- GetAssayData(seurat_split, assay = "RNA", slot = "counts")
  atac_counts_split <- GetAssayData(seurat_split, assay = "ATAC", slot = "counts")
  
  smoothed_split <- smooth_expression(rna_counts_split, atac_counts_split, A_smooth_split)
  cat(sprintf("[OK] %s: %d cells smoothed (preserved all cells)\n", split_name, ncol(smoothed_split$rna_log1p)))
  
  # Save smoothed data
  split_output <- list(
    rna_log1p = smoothed_split$rna_log1p,
    rna_lib = smoothed_split$rna_lib,
    atac_counts = smoothed_split$atac_counts,
    atac_lib = smoothed_split$atac_lib,
    atac_log1p = smoothed_split$atac_log1p,
    rna_scvi_latent = rna_scores_split,
    atac_peakvi_latent = atac_scores_split,
    joint_space = joint_split,
    cells = split_cells,
    split_name = split_name,
    processing_seed = seed,
    k_smooth = k_neighbors,
    method = "scVI_PeakVI"
  )
  saveRDS(split_output, file.path(OUTPUT_METACELLS_DIR, paste0("smoothed_", split_name, ".rds")))
  cat(sprintf("[OK] Saved smoothed_%s.rds\n", split_name))
}

# These Diagnostic Plots can help visualize the results of the smoothing process but they are not included in the manuscript
# ============================================================================
# DIAGNOSTIC PLOTS
# ============================================================================
cat("\n=== Generating diagnostic plots ===\n\n")

# Create plots directory
plots_dir <- file.path(OUTPUT_METACELLS_DIR, "plots_scvi_peakvi")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Load all smoothed data for plotting
smoothed_train <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_train.rds"))
smoothed_val <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_validation.rds"))
smoothed_test <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_test.rds"))

# Use colors specified 
split_colors <- split_colors_pub
modality_colors <- modality_colors_pub

# ----------------------------------------------------------------------------
# PLOT 1: UMAP of Joint Latent Space
# ----------------------------------------------------------------------------
cat("Creating UMAP visualizations...\n")

if (requireNamespace("uwot", quietly = TRUE)) {
  
  # Combine all joint spaces
  all_joint <- rbind(
    smoothed_train$joint_space,
    smoothed_val$joint_space,
    smoothed_test$joint_space
  )
  
  all_splits <- c(
    rep("train", nrow(smoothed_train$joint_space)),
    rep("validation", nrow(smoothed_val$joint_space)),
    rep("test", nrow(smoothed_test$joint_space))
  )
  
  # Run UMAP
  set.seed(seed)
  umap_coords <- uwot::umap(all_joint, n_neighbors = 15, min_dist = 0.3, 
                            n_components = 2, verbose = FALSE)
  
  umap_df <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    Split = factor(all_splits, levels = c("train", "validation", "test"))
  )
  
  p1a <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Split)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = split_colors) +
    labs(title = "UMAP of Joint scVI + PeakVI Space",
         subtitle = sprintf("scVI: %d dims + PeakVI: %d dims = %d total", 
                           train_transforms$rna_latent_dims, 
                           train_transforms$atac_latent_dims,
                           train_transforms$rna_latent_dims + train_transforms$atac_latent_dims),
         x = "UMAP 1", y = "UMAP 2") +
    theme_pub() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  save_plot(p1a, file.path(plots_dir, "01_umap_scvi_peakvi_by_split"), 
                        width = 8, height = 7)
  cat("  ✓ Saved UMAP by split (PNG + PDF)\n")
  
  # Faceted UMAP
  p1b <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Split)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = split_colors) +
    facet_wrap(~Split, ncol = 3) +
    labs(title = "UMAP of scVI + PeakVI Space - Split Comparison",
         x = "UMAP 1", y = "UMAP 2") +
    theme_pub() +
    theme(legend.position = "none")
  
  save  _plot(p1b, file.path(plots_dir, "02_umap_scvi_peakvi_faceted"), 
                        width = 12, height = 5)
  cat("  ✓ Saved UMAP faceted (PNG + PDF)\n")
  
} else {
  cat("  [SKIP] uwot not installed - skipping UMAP plots\n")
}

# ----------------------------------------------------------------------------
# PLOT 2: Latent Variance Distribution (per modality)
# ----------------------------------------------------------------------------
cat("Creating latent variance plots...\n")

# scVI (RNA) variance
rna_var <- apply(smoothed_train$rna_scvi_latent, 2, var)
rna_var_pct <- 100 * rna_var / sum(rna_var)
rna_cumvar <- cumsum(rna_var_pct)

# PeakVI (ATAC) variance
atac_var <- apply(smoothed_train$atac_peakvi_latent, 2, var)
atac_var_pct <- 100 * atac_var / sum(atac_var)
atac_cumvar <- cumsum(atac_var_pct)

var_df <- rbind(
  data.frame(Dim = 1:length(rna_var_pct), Variance = rna_var_pct, 
             Cumulative = rna_cumvar, Modality = "scVI (RNA)"),
  data.frame(Dim = 1:length(atac_var_pct), Variance = atac_var_pct, 
             Cumulative = atac_cumvar, Modality = "PeakVI (ATAC)")
)

p2 <- ggplot(var_df, aes(x = Dim)) +
  geom_bar(aes(y = Variance, fill = Modality), stat = "identity", alpha = 0.8) +
  geom_line(aes(y = Cumulative, color = Modality), linewidth = 1) +
  geom_point(aes(y = Cumulative, color = Modality), size = 2) +
  scale_fill_manual(values = modality_colors) +
  scale_color_manual(values = c("scVI (RNA)" = "#D55E00", "PeakVI (ATAC)" = "#D55E00")) +
  facet_wrap(~Modality, scales = "free_x") +
  scale_y_continuous(
    name = "Variance Explained (%)",
    sec.axis = sec_axis(~., name = "Cumulative Variance (%)")
  ) +
  labs(title = "Latent Dimension Variance by Modality",
       x = "Latent Dimension") +
  theme() +
  theme(legend.position = "none")

save_plot(p2, file.path(plots_dir, "03_latent_variance_by_modality"), 
          width = 10, height = 5)
cat("  ✓ Saved latent variance plot (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 3: Split Cell Counts
# ----------------------------------------------------------------------------
cat("Creating split distribution plots...\n")

split_counts <- data.frame(
  Split = c("train", "validation", "test"),
  Cells = c(length(smoothed_train$cells), 
            length(smoothed_val$cells), 
            length(smoothed_test$cells))
)
split_counts$Percentage <- 100 * split_counts$Cells / sum(split_counts$Cells)
split_counts$Split <- factor(split_counts$Split, levels = c("train", "validation", "test"))

p3 <- ggplot(split_counts, aes(x = Split, y = Cells, fill = Split)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Cells, Percentage)), 
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = split_colors) +
  labs(title = "Cell Distribution Across Data Splits",
       subtitle = sprintf("Total: %d cells", sum(split_counts$Cells)),
       x = NULL, y = "Number of Cells") +
  theme_publication() +
  theme(legend.position = "none") +
  ylim(0, max(split_counts$Cells) * 1.15)

save_plot(p3, file.path(plots_dir, "04_split_distribution"), 
                      width = 6, height = 5)
cat("  ✓ Saved split distribution (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 4: k-NN Distance Distribution
# ----------------------------------------------------------------------------
cat("Creating k-NN distance plots...\n")

nn_train <- RANN::nn2(smoothed_train$joint_space, smoothed_train$joint_space, 
                      k = min(k_neighbors + 1, nrow(smoothed_train$joint_space)))
knn_dist_train <- nn_train$nn.dists[, -1]

nn_val <- RANN::nn2(smoothed_val$joint_space, smoothed_val$joint_space, 
                    k = min(k_neighbors + 1, nrow(smoothed_val$joint_space)))
knn_dist_val <- nn_val$nn.dists[, -1]

nn_test <- RANN::nn2(smoothed_test$joint_space, smoothed_test$joint_space, 
                     k = min(k_neighbors + 1, nrow(smoothed_test$joint_space)))
knn_dist_test <- nn_test$nn.dists[, -1]

dist_df <- data.frame(
  MeanDist = c(rowMeans(knn_dist_train), rowMeans(knn_dist_val), rowMeans(knn_dist_test)),
  Split = c(rep("train", nrow(knn_dist_train)), 
            rep("validation", nrow(knn_dist_val)), 
            rep("test", nrow(knn_dist_test)))
)
dist_df$Split <- factor(dist_df$Split, levels = c("train", "validation", "test"))

p4 <- ggplot(dist_df, aes(x = Split, y = MeanDist, fill = Split)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = split_colors) +
  labs(title = sprintf("k-NN Distance Distribution (k=%d)", k_neighbors),
       subtitle = "Mean distance to k neighbors in scVI+PeakVI space",
       x = NULL, y = "Mean Distance to k Neighbors") +
  theme_pub() +
  theme(legend.position = "none")

save_plot(p4, file.path(plots_dir, "05_knn_distance_distribution"), 
                      width = 7, height = 5)
cat("  ✓ Saved k-NN distance distribution (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 5: Before/After Smoothing - Gene Expression Variance
# ----------------------------------------------------------------------------
cat("Creating smoothing effect plots...\n")

# Get original (unsmoothed) data for training cells
train_cells <- smoothed_train$cells
rna_original_train <- GetAssayData(seurat_obj[, train_cells], assay = "RNA", layer = "counts")
rna_original_train <- log1p(t(t(rna_original_train) / (Matrix::colSums(rna_original_train) + 1e-8)) * 1e6)

# Sample genes for visualization
set.seed(seed)
sample_genes <- sample(rownames(smoothed_train$rna_log1p), 
                       min(1000, nrow(smoothed_train$rna_log1p)))

# Calculate variance per gene
var_original <- apply(rna_original_train[sample_genes, ], 1, var)
var_smoothed <- apply(smoothed_train$rna_log1p[sample_genes, ], 1, var)

var_df2 <- data.frame(
  Original = var_original,
  Smoothed = var_smoothed
)

p5a <- ggplot(var_df2, aes(x = Original, y = Smoothed)) +
  geom_point(alpha = 0.3, size = 1, color = "#0072B2") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Effect of k-NN Smoothing on Gene Expression Variance",
       subtitle = sprintf("Training data (k=%d neighbors, %d genes sampled)", k_neighbors, length(sample_genes)),
       x = "Original Variance (log₁₀)", y = "Smoothed Variance (log₁₀)") +
  theme_pub()

save_plot(p5a, file.path(plots_dir, "06_smoothing_variance_comparison"), 
                      width = 7, height = 6)
cat("  ✓ Saved smoothing variance comparison (PNG + PDF)\n")

# Variance reduction histogram
var_reduction <- (var_original - var_smoothed) / var_original * 100
var_red_df <- data.frame(Reduction = var_reduction[is.finite(var_reduction)])

p5b <- ggplot(var_red_df, aes(x = Reduction)) +
  geom_histogram(bins = 50, fill = "#0072B2", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(var_red_df$Reduction), color = "#D55E00", 
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = median(var_red_df$Reduction), y = Inf, 
           label = sprintf("Median: %.1f%%", median(var_red_df$Reduction)),
           vjust = 2, hjust = -0.1, color = "#D55E00", fontface = "bold", size = 4) +
  labs(title = "Variance Reduction from k-NN Smoothing",
       subtitle = "Positive values indicate variance reduction (smoothing effect)",
       x = "Variance Reduction (%)", y = "Number of Genes") +
  theme_pub()

save_plot(p5b, file.path(plots_dir, "07_variance_reduction_histogram"), 
                      width = 7, height = 5)
cat("  ✓ Saved variance reduction histogram (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 6: Original vs Smoothed - Top Variable Genes
# ----------------------------------------------------------------------------
cat("Creating original vs smoothed gene comparison plots...\n")

top_var_genes <- names(sort(var_original, decreasing = TRUE))[1:6]

scatter_list <- list()
for (i in seq_along(top_var_genes)) {
  gene <- top_var_genes[i]
  gene_df <- data.frame(
    Original = as.numeric(rna_original_train[gene, ]),
    Smoothed = as.numeric(smoothed_train$rna_log1p[gene, ])
  )
  
  scatter_list[[i]] <- ggplot(gene_df, aes(x = Original, y = Smoothed)) +
    geom_point(alpha = 0.3, size = 0.8, color = "#0072B2") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#D55E00", linewidth = 0.5) +
    labs(title = gene, x = "Original", y = "Smoothed") +
    theme_pub(base_size = 10)
}

p6 <- cowplot::plot_grid(plotlist = scatter_list, ncol = 3, nrow = 2)
p6 <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("Original vs Smoothed Expression: Top Variable Genes", 
                                           fontface = "bold", size = 14),
  p6, ncol = 1, rel_heights = c(0.08, 0.92)
)

save_plot(p6, file.path(plots_dir, "08_original_vs_smoothed_top_genes"), 
                      width = 10, height = 7)
cat("  ✓ Saved original vs smoothed top genes (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 7: Expression Density Comparison
# ----------------------------------------------------------------------------
cat("Creating expression density plots...\n")

sorted_var <- sort(var_original, decreasing = TRUE)
gene_high <- names(sorted_var)[1]
gene_med <- names(sorted_var)[floor(length(sorted_var)/2)]
gene_low <- names(sorted_var)[length(sorted_var) - 10]

density_df <- rbind(
  data.frame(Expression = as.numeric(rna_original_train[gene_high, ]), Gene = gene_high, Type = "Original"),
  data.frame(Expression = as.numeric(smoothed_train$rna_log1p[gene_high, ]), Gene = gene_high, Type = "Smoothed"),
  data.frame(Expression = as.numeric(rna_original_train[gene_med, ]), Gene = gene_med, Type = "Original"),
  data.frame(Expression = as.numeric(smoothed_train$rna_log1p[gene_med, ]), Gene = gene_med, Type = "Smoothed"),
  data.frame(Expression = as.numeric(rna_original_train[gene_low, ]), Gene = gene_low, Type = "Original"),
  data.frame(Expression = as.numeric(smoothed_train$rna_log1p[gene_low, ]), Gene = gene_low, Type = "Smoothed")
)

density_df$Gene <- factor(density_df$Gene, levels = c(gene_high, gene_med, gene_low),
                          labels = c(paste0(gene_high, "\n(High var)"),
                                     paste0(gene_med, "\n(Med var)"),
                                     paste0(gene_low, "\n(Low var)")))

p7 <- ggplot(density_df, aes(x = Expression, fill = Type, color = Type)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  facet_wrap(~Gene, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Original" = "#D55E00", "Smoothed" = "#0072B2")) +
  scale_color_manual(values = c("Original" = "#D55E00", "Smoothed" = "#0072B2")) +
  labs(title = "Expression Density: Original vs Smoothed",
       subtitle = "Smoothing reduces sparsity and noise",
       x = "log₁p(CPM)", y = "Density") +
  theme_pub() +
  theme(legend.position = "bottom")

save_plot(p7, file.path(plots_dir, "09_expression_density_comparison"), 
                      width = 11, height = 5)
cat("  ✓ Saved expression density comparison (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 8: Dropout Reduction
# ----------------------------------------------------------------------------
cat("Creating dropout reduction plot...\n")

zero_prop_original <- rowMeans(rna_original_train[sample_genes, ] == 0)
zero_prop_smoothed <- rowMeans(smoothed_train$rna_log1p[sample_genes, ] == 0)

zero_df <- data.frame(
  Original = zero_prop_original,
  Smoothed = zero_prop_smoothed
)

p8 <- ggplot(zero_df, aes(x = Original, y = Smoothed)) +
  geom_point(alpha = 0.3, size = 1, color = "#009E73") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  labs(title = "Dropout Reduction from k-NN Smoothing",
       subtitle = "Proportion of zero-expression cells per gene",
       x = "Original Zero Proportion", y = "Smoothed Zero Proportion") +
  theme_pub() +
  annotate("text", x = 0.8, y = 0.2, 
           label = sprintf("Mean reduction:\n%.1f%% → %.1f%%", 
                           mean(zero_prop_original)*100, mean(zero_prop_smoothed)*100),
           hjust = 0.5, size = 4, fontface = "bold", color = "#009E73")

save_plot(p8, file.path(plots_dir, "10_dropout_reduction"), 
                      width = 7, height = 6)
cat("  ✓ Saved dropout reduction plot (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 9: Mean Expression Preservation
# ----------------------------------------------------------------------------
cat("Creating mean expression preservation plot...\n")

mean_original <- rowMeans(rna_original_train[sample_genes, ])
mean_smoothed <- rowMeans(smoothed_train$rna_log1p[sample_genes, ])

mean_df <- data.frame(
  Original = mean_original,
  Smoothed = mean_smoothed
)

cor_mean <- cor(mean_original, mean_smoothed)

p9 <- ggplot(mean_df, aes(x = Original, y = Smoothed)) +
  geom_point(alpha = 0.3, size = 1, color = "#E69F00") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#0072B2", linewidth = 0.8) +
  labs(title = "Mean Expression Preservation",
       subtitle = sprintf("Correlation: r = %.4f (smoothing preserves mean structure)", cor_mean),
       x = "Original Mean Expression", y = "Smoothed Mean Expression") +
  theme_pub()

save_plot(p9, file.path(plots_dir, "11_mean_expression_preservation"), 
                      width = 7, height = 6)
cat("  ✓ Saved mean expression preservation plot (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 10: Expression Heatmap (Original vs Smoothed)
# ----------------------------------------------------------------------------
cat("Creating expression heatmaps...\n")

top50_genes <- names(sort(var_original, decreasing = TRUE))[1:50]
sample_cells_hm <- sample(train_cells, min(200, length(train_cells)))

orig_mat <- as.matrix(rna_original_train[top50_genes, sample_cells_hm])
smooth_mat <- as.matrix(smoothed_train$rna_log1p[top50_genes, sample_cells_hm])

orig_scaled <- t(scale(t(orig_mat)))
smooth_scaled <- t(scale(t(smooth_mat)))

orig_hm_df <- reshape2::melt(orig_scaled)
colnames(orig_hm_df) <- c("Gene", "Cell", "Expression")
orig_hm_df$Type <- "Original"

smooth_hm_df <- reshape2::melt(smooth_scaled)
colnames(smooth_hm_df) <- c("Gene", "Cell", "Expression")
smooth_hm_df$Type <- "Smoothed"

hm_df <- rbind(orig_hm_df, smooth_hm_df)
hm_df$Type <- factor(hm_df$Type, levels = c("Original", "Smoothed"))
hm_df$Expression <- pmax(pmin(hm_df$Expression, 3), -3)

p10 <- ggplot(hm_df, aes(x = Cell, y = Gene, fill = Expression)) +
  geom_tile() +
  facet_wrap(~Type, ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, limits = c(-3, 3),
                       name = "Scaled\nExpression") +
  labs(title = "Expression Heatmap: Original vs Smoothed",
       subtitle = sprintf("Top 50 variable genes × %d cells (z-scored)", length(sample_cells_hm)),
       x = "Cells", y = "Genes") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 5))

save_plot(p10, file.path(plots_dir, "12_expression_heatmap_comparison"), 
                      width = 14, height = 10)
cat("  ✓ Saved expression heatmap comparison (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 11: scVI-PeakVI Latent Correlation
# ----------------------------------------------------------------------------
cat("Creating scVI-PeakVI correlation plots...\n")

n_comp_plot <- min(5, train_transforms$rna_latent_dims, train_transforms$atac_latent_dims)

rna_latent_train <- smoothed_train$rna_scvi_latent[, 1:n_comp_plot]
atac_latent_train <- smoothed_train$atac_peakvi_latent[, 1:n_comp_plot]

cor_matrix <- cor(rna_latent_train, atac_latent_train)

cor_df <- expand.grid(
  scVI = paste0("scVI_", 1:n_comp_plot),
  PeakVI = paste0("PeakVI_", 1:n_comp_plot)
)
cor_df$Correlation <- as.vector(cor_matrix)

p11 <- ggplot(cor_df, aes(x = scVI, y = PeakVI, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, limits = c(-1, 1),
                       name = "Correlation") +
  labs(title = "scVI vs PeakVI Latent Correlation",
       subtitle = "Cross-modality relationships in joint space",
       x = "scVI (RNA) Dimension", y = "PeakVI (ATAC) Dimension") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p11, file.path(plots_dir, "13_scvi_peakvi_correlation_heatmap"), 
                      width = 7, height = 6)
cat("  ✓ Saved scVI-PeakVI correlation heatmap (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 12: Summary Parameters Panel
# ----------------------------------------------------------------------------
cat("Creating summary panel...\n")

summary_text <- paste0(
  "METACELL CREATION PARAMETERS\n",
  "Method: scVI + PeakVI Latent Space\n\n",
  "Latent Space:\n",
  sprintf("  - scVI (RNA) latent dims: %d\n", train_transforms$rna_latent_dims),
  sprintf("  - PeakVI (ATAC) latent dims: %d\n", train_transforms$atac_latent_dims),
  sprintf("  - Total joint dims: %d\n", train_transforms$rna_latent_dims + train_transforms$atac_latent_dims),
  sprintf("  - Z-score normalized: Yes\n\n"),
  "k-NN Smoothing:\n",
  sprintf("  - k neighbors: %d\n", k_neighbors),
  sprintf("  - Random seed: %d\n\n", seed),
  "Data Splits:\n",
  sprintf("  - Training cells: %d\n", length(smoothed_train$cells)),
  sprintf("  - Validation cells: %d\n", length(smoothed_val$cells)),
  sprintf("  - Test cells: %d\n", length(smoothed_test$cells)),
  sprintf("  - Total cells: %d", 
          length(smoothed_train$cells) + length(smoothed_val$cells) + length(smoothed_test$cells))
)

p12 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = summary_text, 
           hjust = 0.5, vjust = 0.5, size = 4, family = "mono") +
  theme_void() +
  theme(plot.margin = margin(20, 20, 20, 20),
        plot.background = element_rect(fill = "white", color = NA))

save_plot(p12, file.path(plots_dir, "14_summary_parameters"), 
                      width = 6, height = 5)
cat("  ✓ Saved summary parameters (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# COMBINED PDF REPORT
# ----------------------------------------------------------------------------
cat("\nCreating combined PDF report...\n")

cairo_pdf(file.path(plots_dir, "Step_050_Metacell_scVI_PeakVI_Report.pdf"), width = 10, height = 8)

# Title page
plot.new()
text(0.5, 0.6, "Step 050: Metacell Creation Report", cex = 2, font = 2)
text(0.5, 0.5, "scVI + PeakVI Method", cex = 1.6, font = 3)
text(0.5, 0.4, sprintf("Sample: %s", SAMPLE_NAME), cex = 1.4)
text(0.5, 0.3, sprintf("Date: %s", Sys.Date()), cex = 1.2)
text(0.5, 0.2, "Figures", cex = 1, font = 3)

# Print all plots
if (exists("p1a")) print(p1a)
if (exists("p1b")) print(p1b)
print(p2)
print(p3)
print(p4)
print(p5a)
print(p5b)
print(p6)
print(p7)
print(p8)
print(p9)
print(p10)
print(p11)
print(p12)

dev.off()
cat("  ✓ Saved combined PDF report (vector graphics)\n")

cat("\n")
cat("plots saved to:", plots_dir, "\n")
cat("  - Each plot saved as both .png (600 dpi) and .pdf (vector)\n")
cat("  - Combined report: Step_030_Metacell_scVI_PeakVI_Report.pdf\n")

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 050 COMPLETE (scVI + PeakVI)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("k-NN smoothing completed for all data splits using scVI + PeakVI latent space\n\n")

cat("Parameters:\n")
cat(sprintf("  - k = %d neighbors\n", k_neighbors))
cat(sprintf("  - scVI RNA latent dims = %d\n", train_transforms$rna_latent_dims))
cat(sprintf("  - PeakVI ATAC latent dims = %d\n", train_transforms$atac_latent_dims))
cat(sprintf("  - Total joint dims = %d\n", train_transforms$rna_latent_dims + train_transforms$atac_latent_dims))
cat(sprintf("  - Random seed = %d\n\n", seed))

cat("Key features:\n")
cat("  ✓ Training: Used scVI (RNA) + PeakVI (ATAC) concatenated latent space\n")
cat("  ✓ Val/Test: Applied training-based z-score to each modality's latent\n")
cat("  ✓ Val/Test: Smoothed only within each split (no leakage)\n")
cat("  ✓ All cells preserved (1-to-1 mapping, no aggregation)\n")
cat("  ✓ plots (PNG + PDF)\n\n")

cat("Output files:\n")
cat("  - train_transforms_scvi_peakvi.rds (z-score parameters)\n")
cat("  - smoothed_train.rds\n")
cat("  - smoothed_validation.rds\n")
cat("  - smoothed_test.rds\n\n")

cat("Diagnostic plots (PNG + PDF formats):\n")
cat("  01 - UMAP by split\n")
cat("  02 - UMAP faceted\n")
cat("  03 - Latent variance by modality\n")
cat("  04 - Split distribution\n")
cat("  05 - k-NN distance distribution\n")
cat("  06 - Smoothing variance comparison\n")
cat("  07 - Variance reduction histogram\n")
cat("  08 - Original vs smoothed top genes\n")
cat("  09 - Expression density comparison\n")
cat("  10 - Dropout reduction\n")
cat("  11 - Mean expression preservation\n")
cat("  12 - Expression heatmap comparison\n")
cat("  13 - scVI-PeakVI correlation heatmap\n")
cat("  14 - Summary parameters\n")
cat("  + Combined PDF report\n\n")

cat("Output directory:", OUTPUT_METACELLS_DIR, "\n")
cat("Plots directory:", plots_dir, "\n\n")

cat("Next step: Run 04_feature_extraction.R\n")
