# ============================================================================
# Step_030: Metacell Creation via KNN Smoothing (Automated)
# ============================================================================
# This script creates smoothed expression profiles using k-NN smoothing
# in a joint RNA+ATAC latent space. It reads configuration from config.R.
#
# APPROACH:
# - Project RNA and ATAC into low-dimensional spaces (PCA and LSI)
# - Combine into joint multimodal space with z-score standardization
# - For each cell, find k nearest neighbors in joint space
# - Smooth expression by averaging with k neighbors (preserves all cells)
# - Training: build PCA/LSI transformations
# - Val/Test: project into training space, smooth among val/test cells only
#
# Input: Seurat object with split labels from Step_020
# Output: Smoothed expression matrices for train/val/test splits
#
# Usage:
#   1. Ensure config.R is properly configured
#   2. Run: Rscript 03_metacell_creation_pca_lsi.R
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
  config_path <- "config.R"
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
  library(Signac)
  library(Matrix)
  library(ggplot2)
  library(RANN)
  library(dplyr)
  library(cowplot)
  library(reshape2)
})

# Check for irlba
if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("Package 'irlba' is required. Install with: install.packages('irlba')")
}

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

#' Plotting ggplot2 theme
#' @param base_size Base font size (default 12)
#' @param base_family Font family (default "Arial")
theme_pub <- function(base_size = 12, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Text elements
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5, 
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "grey30",
                                   margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 2, hjust = 1, color = "grey50"),
      
      # Axis elements
      axis.title = element_text(face = "bold", size = base_size),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      
      # Legend elements
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.background = element_blank(),
      legend.key = element_blank(),
      
      # Panel elements
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Facet elements
      strip.text = element_text(face = "bold", size = base_size, color = "black"),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      
      # Plot margins
      plot.margin = margin(15, 15, 15, 15),
      
      # Plot tag (for multi-panel figures)
      plot.tag = element_text(face = "bold", size = base_size + 4)
    )
}

#' Colorblind-friendly palette
colorblind_palette <- c(
  "#E69F00",  # Orange

"#56B4E9",  # Sky Blue
  "#009E73",  # Bluish Green
  "#F0E442",  # Yellow
  "#0072B2",  # Blue
  "#D55E00",  # Vermillion
  "#CC79A7",  # Reddish Purple
  "#999999"   # Grey
)

#' Split-specific colors (consistent across all plots)
split_colors_pub <- c(
  "train" = "#0072B2",       # Blue
  "validation" = "#E69F00",  # Orange  
  "test" = "#009E73"         # Green
)

#' Save plot in both PNG and PDF formats
#' @param plot ggplot object
#' @param filename Base filename without extension
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution for PNG
save_plot <- function(plot, filename, width = 7, height = 5, dpi = 600) {
  # Save PNG 
  ggsave(paste0(filename, ".png"), plot = plot, 
         width = width, height = height, dpi = dpi)
  
 # Save PDF 
  ggsave(paste0(filename, ".pdf"), plot = plot, 
         width = width, height = height, device = cairo_pdf)
  
  invisible(plot)
}

# Create output directories
create_output_directories()

# Set parameters from config
pca_dims <- PCA_DIMS
k_neighbors <- K_NEIGHBORS
n_variable_features <- N_VARIABLE_FEATURES
seed <- SEED_METACELL

set.seed(seed)

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 030: Metacell Creation (k-NN Smoothing)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Parameters:\n")
cat(sprintf("  PCA/LSI dimensions: %d\n", pca_dims))
cat(sprintf("  k-NN neighbors: %d\n", k_neighbors))
cat(sprintf("  Variable features: %d\n", n_variable_features))
cat(sprintf("  Random seed: %d\n\n", seed))

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

# Normalize RNA counts
normalize_rna <- function(obj) { 
  DefaultAssay(obj) <- "RNA"
  NormalizeData(obj, verbose = FALSE)
}

# Build PCA transformation from training data
build_rna_pca_fit <- function(obj, pca_dims, nfeatures = 3000) {
  obj <- normalize_rna(obj)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE)
  vf  <- VariableFeatures(obj)
  
  X   <- GetAssayData(obj, assay = "RNA", layer = "data")[vf, , drop = FALSE]
  mu  <- Matrix::rowMeans(X)
  Xsq <- X
  if (inherits(Xsq, "dgCMatrix")) Xsq@x <- Xsq@x^2 else Xsq <- X^2
  v   <- Matrix::rowMeans(Xsq) - mu^2
  sd  <- sqrt(pmax(1e-8, v))
  
  Z   <- t((X - mu) / sd)
  pr  <- irlba::prcomp_irlba(Z, n = pca_dims, center = FALSE, scale. = FALSE)
  
  R <- pr$rotation
  if (is.null(rownames(R))) rownames(R) <- vf
  feat_names <- rownames(R)
  
  list(
    vf = vf,
    mu = as.numeric(mu), names_mu = names(mu),
    sd = as.numeric(sd), names_sd = names(sd),
    rotation = R,
    feat_names = feat_names,
    train_scores = pr$x
  )
}

# Project new data into training PCA space
project_rna_pca <- function(obj, fit) {
  X_all <- GetAssayData(obj, assay = "RNA", layer = "data")
  
  feat_names <- if (!is.null(fit$feat_names)) fit$feat_names
  else if (!is.null(rownames(fit$rotation))) rownames(fit$rotation)
  else fit$vf
  
  common <- intersect(feat_names, rownames(X_all))
  if (length(common) < max(50, floor(0.5 * length(feat_names)))) {
    stop(sprintf("Too few overlapping genes: %d of %d", length(common), length(feat_names)))
  }
  
  ord <- match(common, feat_names)
  common <- common[order(ord)]
  
  Rot <- fit$rotation
  if (is.null(rownames(Rot))) rownames(Rot) <- feat_names
  
  X <- X_all[common, , drop = FALSE]
  R <- Rot[common, , drop = FALSE]
  
  mu <- fit$mu; names(mu) <- fit$names_mu
  sd <- fit$sd; names(sd) <- fit$names_sd
  mu <- mu[common]; sd <- sd[common]
  sd[!is.finite(sd) | sd == 0] <- 1
  
  Z <- t((X - mu) / sd)
  as.matrix(Z %*% R)
}

# Build LSI transformation from TRAINING ATAC data
build_atac_lsi_signac <- function(seurat_obj, assay = "ATAC", ncomp = 30) {
  DefaultAssay(seurat_obj) <- assay
  
  counts <- GetAssayData(seurat_obj, assay = assay, slot = "counts")
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
    seurat_obj <- SetAssayData(seurat_obj, assay = assay, slot = "counts", new.data = counts)
  }
  
  cat("[INFO] Selecting top features for ATAC data (min.cutoff = 'q5')...\n")
  seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = "q05")
  
  cat("[INFO] Running Signac TF-IDF normalization...\n")
  seurat_obj <- RunTFIDF(seurat_obj, verbose = FALSE)
  
  n_cells_train <- ncol(counts)
  nz_per_peak   <- Matrix::rowSums(counts > 0)
  idf_vector    <- log1p(n_cells_train / pmax(1, nz_per_peak))
  names(idf_vector) <- rownames(counts)
  
  cat(sprintf("[INFO] Running Signac SVD (LSI) with %d components...\n", ncomp))
  seurat_obj <- RunSVD(seurat_obj, n = ncomp, verbose = FALSE)
  
  lsi_coords_full       <- Embeddings(seurat_obj, "lsi")
  feature_loadings_full <- Loadings(seurat_obj, "lsi")
  
  lsi_coords      <- lsi_coords_full[, -1, drop = FALSE]
  feature_loadings <- feature_loadings_full[, -1, drop = FALSE]
  
  cat(sprintf("[INFO] ATAC LSI complete: %d peaks × %d cells; %d components (skipped 1st)\n",
              nrow(feature_loadings), nrow(lsi_coords), ncol(lsi_coords)))
  
  list(
    seurat_obj       = seurat_obj,
    feature_loadings = feature_loadings,
    train_scores     = lsi_coords,
    n_components     = ncol(lsi_coords),
    idf              = idf_vector
  )
}

# Project VAL/TEST (new data) split into TRAIN LSI space
project_atac_lsi_signac <- function(seurat_obj_split, train_fit, assay = "ATAC") {
  DefaultAssay(seurat_obj_split) <- assay
  
  cat("[INFO] Projecting ATAC data into TRAINING LSI space...\n")
  
  L      <- train_fit$feature_loadings
  idf    <- train_fit$idf
  tpeaks <- rownames(L)
  
  Xc <- GetAssayData(seurat_obj_split, assay = assay, slot = "counts")
  if (!inherits(Xc, "dgCMatrix")) {
    Xc <- as(Xc, "dgCMatrix")
  }
  
  common <- tpeaks[tpeaks %in% rownames(Xc)]
  if (length(common) < max(500, floor(0.5 * length(tpeaks)))) {
    stop(sprintf("Too few overlapping peaks for projection: %d/%d", length(common), length(tpeaks)))
  }
  
  Xc   <- Xc[common, , drop = FALSE]
  Luse <- L[common, , drop = FALSE]
  idf  <- idf[common]
  
  lib <- Matrix::colSums(Xc)
  lib[lib == 0] <- 1
  TF <- t(t(Xc) / lib)
  
  TFIDF <- TF * idf
  
  emb <- as.matrix(Matrix::t(TFIDF) %*% Luse)
  rownames(emb) <- colnames(seurat_obj_split)
  
  cat(sprintf("[INFO] Projected %d cells into %d-dimensional ATAC LSI space (excl. comp 1)\n",
              nrow(emb), ncol(emb)))
  
  emb
}

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

# Build RNA PCA from training data
cat("\n[INFO] Building RNA PCA from training data...\n")
rna_fit <- build_rna_pca_fit(seurat_train, pca_dims = pca_dims, nfeatures = n_variable_features)
rna_scores_train <- rna_fit$train_scores
rownames(rna_scores_train) <- train_cells

cat(sprintf("RNA PCA: %d variable genes, %d components\n", length(rna_fit$vf), ncol(rna_scores_train)))

# Build ATAC LSI from training data
cat("\n[INFO] Building ATAC LSI from training data...\n")
atac_fit <- build_atac_lsi_signac(seurat_train, assay = "ATAC", ncomp = pca_dims)
atac_scores_train <- atac_fit$train_scores
rownames(atac_scores_train) <- train_cells

cat(sprintf("ATAC LSI: %d peaks, %d components\n", nrow(atac_fit$feature_loadings), ncol(atac_scores_train)))

# Z-score standardization
cat("\n[INFO] Computing z-score parameters from training...\n")
rna_zs <- z_fit_cols(rna_scores_train)
atac_zs <- z_fit_cols(atac_scores_train)

rna_train_z <- z_apply_cols(rna_scores_train, rna_zs$mu, rna_zs$sd)
atac_train_z <- z_apply_cols(atac_scores_train, atac_zs$mu, atac_zs$sd)

# Create joint space
joint_train <- cbind(rna_train_z, atac_train_z)
cat(sprintf("Joint space: %d cells × %d dimensions\n", nrow(joint_train), ncol(joint_train)))

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
  rna_fit = rna_fit,
  atac_fit = atac_fit,
  rna_zs = rna_zs,
  atac_zs = atac_zs,
  pca_dims = pca_dims,
  k_neighbors = k_neighbors,
  n_variable_features = n_variable_features
)
saveRDS(train_transforms, file.path(OUTPUT_METACELLS_DIR, "train_transforms.rds"))
cat("[OK] Saved train_transforms.rds\n")

# Save smoothed training data
train_output <- list(
  rna_log1p = smoothed_train$rna_log1p,
  rna_lib = smoothed_train$rna_lib,
  atac_counts = smoothed_train$atac_counts,
  atac_lib = smoothed_train$atac_lib,
  atac_log1p = smoothed_train$atac_log1p,
  rna_pca = rna_scores_train,
  atac_lsi = atac_scores_train,
  joint_space = joint_train,
  cells = train_cells,
  split_name = "train",
  processing_seed = seed,
  k_smooth = k_neighbors,
  pca_dims = pca_dims
)
saveRDS(train_output, file.path(OUTPUT_METACELLS_DIR, "smoothed_train.rds"))
cat("[OK] Saved smoothed_train.rds\n")

# ============================================================================
# PROCESS VALIDATION AND TEST DATA
# ============================================================================
cat("\n=== Processing VALIDATION and TEST data ===\n")

train_transforms <- readRDS(file.path(OUTPUT_METACELLS_DIR, "train_transforms.rds"))

for (split_name in c("validation", "test")) {
  cat(sprintf("\n--- Processing %s ---\n", toupper(split_name)))
  
  split_cells <- colnames(seurat_obj)[seurat_obj$data_split == split_name]
  seurat_split <- seurat_obj[, split_cells]
  
  cat(sprintf("%s cells: %d\n", split_name, length(split_cells)))
  
  # Ensure data is sparse and normalized
  DefaultAssay(seurat_split) <- "RNA"
  seurat_split <- normalize_rna(seurat_split)
  
  rc <- GetAssayData(seurat_split, assay = "RNA", layer = "counts")
  if (!inherits(rc, "dgCMatrix")) {
    seurat_split <- SetAssayData(seurat_split, assay = "RNA", layer = "counts", new.data = as(rc, "dgCMatrix"))
  }
  
  ac <- GetAssayData(seurat_split, assay = "ATAC", layer = "counts")
  if (!inherits(ac, "dgCMatrix")) {
    seurat_split <- SetAssayData(seurat_split, assay = "ATAC", layer = "counts", new.data = as(ac, "dgCMatrix"))
  }
  
  # PROJECT into training PCA space
  cat(sprintf("[INFO] Projecting %s RNA into TRAINING PCA space...\n", split_name))
  rna_scores_split <- project_rna_pca(seurat_split, train_transforms$rna_fit)
  rownames(rna_scores_split) <- split_cells
  cat(sprintf("  Projected using %d genes\n", nrow(train_transforms$rna_fit$rotation)))
  
  # PROJECT into training LSI space
  cat(sprintf("[INFO] Projecting %s ATAC into TRAINING LSI space...\n", split_name))
  atac_scores_split <- project_atac_lsi_signac(seurat_split, train_transforms$atac_fit, assay = "ATAC")
  rownames(atac_scores_split) <- split_cells
  cat(sprintf("  Projected using %d peaks\n", nrow(train_transforms$atac_fit$feature_loadings)))
  
  # Apply TRAINING z-score parameters
  cat("[INFO] Applying training z-score normalization...\n")
  rna_split_z <- z_apply_cols(rna_scores_split, train_transforms$rna_zs$mu, train_transforms$rna_zs$sd)
  atac_split_z <- z_apply_cols(atac_scores_split, train_transforms$atac_zs$mu, train_transforms$atac_zs$sd)
  
  # Create joint space
  joint_split <- cbind(rna_split_z, atac_split_z)
  cat(sprintf("Joint space: %d cells × %d dimensions\n", nrow(joint_split), ncol(joint_split)))
  
  # Build k-NN smoothing matrix WITHIN this split (prevents data leakage)
  cat(sprintf("[INFO] Building k-NN within %s cells (k=%d, no training data)...\n", split_name, k_neighbors))
  A_smooth_split <- build_knn_smooth(joint_split, split_cells, k_smooth = k_neighbors, seed = seed)
  
  neighbors_per_cell <- Matrix::rowSums(A_smooth_split > 0)
  cat(sprintf("k-NN stats: mean=%.1f, min=%d, max=%d neighbors\n",
              mean(neighbors_per_cell), min(neighbors_per_cell), max(neighbors_per_cell)))
  
  # Smooth using ONLY cells from this split
  cat(sprintf("[INFO] Smoothing %s expression...\n", split_name))
  rna_counts_split <- GetAssayData(seurat_split, assay = "RNA", layer = "counts")
  atac_counts_split <- GetAssayData(seurat_split, assay = "ATAC", layer = "counts")
  
  smoothed_split <- smooth_expression(rna_counts_split, atac_counts_split, A_smooth_split)
  cat(sprintf("[OK] %s: %d cells smoothed (preserved all cells)\n", split_name, ncol(smoothed_split$rna_log1p)))
  
  # Save smoothed (metacell) data
  split_output <- list(
    rna_log1p = smoothed_split$rna_log1p,
    rna_lib = smoothed_split$rna_lib,
    atac_counts = smoothed_split$atac_counts,
    atac_lib = smoothed_split$atac_lib,
    atac_log1p = smoothed_split$atac_log1p,
    rna_pca = rna_scores_split,
    atac_lsi = atac_scores_split,
    joint_space = joint_split,
    cells = split_cells,
    split_name = split_name,
    processing_seed = seed,
    k_smooth = k_neighbors,
    pca_dims = pca_dims
  )
  saveRDS(split_output, file.path(OUTPUT_METACELLS_DIR, paste0("smoothed_", split_name, ".rds")))
  cat(sprintf("[OK] Saved smoothed_%s.rds\n", split_name))
}

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================
# These Diagnostic plots help visualize the metacell creation process. They are not included in the manuscript but provide insights into the data structure and smoothing effects.
# ============================================================================
# DIAGNOSTIC PLOTS 
# ============================================================================
cat("\n=== Generating diagnostic plots ===\n\n")

# Create plots directory
plots_dir <- file.path(OUTPUT_METACELLS_DIR, "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Load all smoothed (metacell) data for plotting
smoothed_train <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_train.rds"))
smoothed_val <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_validation.rds"))
smoothed_test <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_test.rds"))

# Use publication-ready split colors (defined in theme section)
split_colors <- split_colors_pub

# ----------------------------------------------------------------------------
# PLOT 1: UMAP of Joint Space (if uwot is available)
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
  
  # Run UMAP on joint space
  set.seed(seed)
  umap_coords <- uwot::umap(all_joint, n_neighbors = 15, min_dist = 0.3, 
                            n_components = 2, verbose = FALSE)
  
  umap_df <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    Split = factor(all_splits, levels = c("train", "validation", "test"))
  )
  
  # Plot 1a: UMAP colored by split
  p1a <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Split)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = split_colors) +
    labs(title = "UMAP of Joint RNA+ATAC Space",
         subtitle = sprintf("PCA: %d dims, LSI: %d dims (z-scored)", pca_dims, pca_dims),
         x = "UMAP 1", y = "UMAP 2") +
    theme_pub() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  save_plot(p1a, file.path(plots_dir, "01_umap_joint_space_by_split"), 
                        width = 8, height = 7)
  cat("  ✓ Saved UMAP by split (PNG + PDF)\n")
  
  # Plot 1b: UMAP split into facets
  p1b <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Split)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = split_colors) +
    facet_wrap(~Split, ncol = 3) +
    labs(title = "UMAP of Joint Space - Split Comparison",
         x = "UMAP 1", y = "UMAP 2") +
    theme_pub() +
    theme(legend.position = "none")
  
  save_plot(p1b, file.path(plots_dir, "02_umap_joint_space_faceted"), 
                        width = 12, height = 5)
  cat("  ✓ Saved UMAP faceted\n")
  
} else {
  cat("  [SKIP] uwot not installed - skipping UMAP plots\n")
}

# ----------------------------------------------------------------------------
# PLOT 2: PCA Variance Explained
# ----------------------------------------------------------------------------
cat("Creating PCA/LSI variance plots...\n")

# Get variance from PCA (if available in fit)
if (!is.null(train_transforms$rna_fit$train_scores)) {
  pca_var <- apply(train_transforms$rna_fit$train_scores, 2, var)
  pca_var_pct <- 100 * pca_var / sum(pca_var)
  pca_cumvar <- cumsum(pca_var_pct)
  
  pca_df <- data.frame(
    PC = 1:length(pca_var_pct),
    Variance = pca_var_pct,
    Cumulative = pca_cumvar
  )
  
  p2a <- ggplot(pca_df, aes(x = PC)) +
    geom_bar(aes(y = Variance), stat = "identity", fill = "#56B4E9", alpha = 0.8) +
    geom_line(aes(y = Cumulative), color = "#D55E00", linewidth = 1) +
    geom_point(aes(y = Cumulative), color = "#D55E00", size = 2) +
    scale_y_continuous(
      name = "Variance Explained (%)",
      sec.axis = sec_axis(~., name = "Cumulative Variance (%)")
    ) +
    labs(title = "RNA PCA Variance Explained",
         subtitle = sprintf("First %d components", pca_dims),
         x = "Principal Component") +
    theme_pub()
  
  save_plot(p2a, file.path(plots_dir, "03_pca_variance_explained"), 
                        width = 8, height = 5)
  cat("  ✓ Saved PCA variance plot (PNG + PDF)\n")
}

# LSI variance
if (!is.null(train_transforms$atac_fit$train_scores)) {
  lsi_var <- apply(train_transforms$atac_fit$train_scores, 2, var)
  lsi_var_pct <- 100 * lsi_var / sum(lsi_var)
  lsi_cumvar <- cumsum(lsi_var_pct)
  
  lsi_df <- data.frame(
    LSI = 1:length(lsi_var_pct),
    Variance = lsi_var_pct,
    Cumulative = lsi_cumvar
  )
  
  p2b <- ggplot(lsi_df, aes(x = LSI)) +
    geom_bar(aes(y = Variance), stat = "identity", fill = "#009E73", alpha = 0.8) +
    geom_line(aes(y = Cumulative), color = "#D55E00", linewidth = 1) +
    geom_point(aes(y = Cumulative), color = "#D55E00", size = 2) +
    scale_y_continuous(
      name = "Variance Explained (%)",
      sec.axis = sec_axis(~., name = "Cumulative Variance (%)")
    ) +
    labs(title = "ATAC LSI Variance Explained",
         subtitle = sprintf("Components 2-%d (1st excluded)", pca_dims),
         x = "LSI Component") +
    theme_pub()
  
  save_plot(p2b, file.path(plots_dir, "04_lsi_variance_explained"), 
                        width = 8, height = 5)
  cat("  ✓ Saved LSI variance plot (PNG + PDF)\n")
}

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
  theme_pub() +
  theme(legend.position = "none") +
  ylim(0, max(split_counts$Cells) * 1.15)

save_plot(p3, file.path(plots_dir, "05_split_distribution"), 
          width = 6, height = 5)
cat("  ✓ Saved split distribution (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 4: k-NN Distance Distribution
# ----------------------------------------------------------------------------
cat("Creating k-NN distance plots...\n")

# Compute k-NN distances for training data
nn_train <- RANN::nn2(smoothed_train$joint_space, smoothed_train$joint_space, 
                      k = min(k_neighbors + 1, nrow(smoothed_train$joint_space)))
knn_dist_train <- nn_train$nn.dists[, -1]  # Exclude self

# Sample for validation/test
nn_val <- RANN::nn2(smoothed_val$joint_space, smoothed_val$joint_space, 
                    k = min(k_neighbors + 1, nrow(smoothed_val$joint_space)))
knn_dist_val <- nn_val$nn.dists[, -1]

nn_test <- RANN::nn2(smoothed_test$joint_space, smoothed_test$joint_space, 
                     k = min(k_neighbors + 1, nrow(smoothed_test$joint_space)))
knn_dist_test <- nn_test$nn.dists[, -1]

# Create distance dataframe (mean distance to k neighbors per cell)
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
       subtitle = "Mean Euclidean distance to k nearest neighbors in joint space",
       x = NULL, y = "Mean Distance to k Neighbors") +
  theme_pub() +
  theme(legend.position = "none")

save_plot(p4, file.path(plots_dir, "06_knn_distance_distribution"), 
                      width = 7, height = 5)
cat("  ✓ Saved k-NN distance distribution (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 5: Before/After Smoothing - Gene Expression Variance
# ----------------------------------------------------------------------------
cat("Creating smoothing effect plots...\n")

# Sample some genes for visualization
set.seed(seed)
sample_genes <- sample(rownames(smoothed_train$rna_log1p), min(1000, nrow(smoothed_train$rna_log1p)))

# Get original (unsmoothed) data for comparison
rna_original_train <- GetAssayData(seurat_obj[, smoothed_train$cells], assay = "RNA", layer = "counts")
rna_original_train <- log1p(t(t(rna_original_train) / (Matrix::colSums(rna_original_train) + 1e-8)) * 1e6)

# Calculate variance per gene
var_original <- apply(rna_original_train[sample_genes, ], 1, var)
var_smoothed <- apply(smoothed_train$rna_log1p[sample_genes, ], 1, var)

var_df <- data.frame(
  Original = var_original,
  Smoothed = var_smoothed
)

p5a <- ggplot(var_df, aes(x = Original, y = Smoothed)) +
  geom_point(alpha = 0.3, size = 1, color = "#0072B2") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Effect of k-NN Smoothing on Gene Expression Variance",
       subtitle = sprintf("Training data (k=%d neighbors, %d genes sampled)", k_neighbors, length(sample_genes)),
       x = "Original Variance (log₁₀)", y = "Smoothed Variance (log₁₀)") +
  theme_pub()

save_plot(p5a, file.path(plots_dir, "07_smoothing_variance_comparison"), 
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

save_plot(p5b, file.path(plots_dir, "08_variance_reduction_histogram"), 
                      width = 7, height = 5)
cat("  ✓ Saved variance reduction histogram (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 5c-5g: Additional Original vs Smoothed Comparisons
# ----------------------------------------------------------------------------
cat("Creating additional original vs smoothed comparison plots...\n")

# PLOT 5c: Side-by-side scatter for top variable genes (individual gene examples)
# Select top 6 most variable genes for visualization
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

# Combine into grid
p5c <- cowplot::plot_grid(plotlist = scatter_list, ncol = 3, nrow = 2)
p5c <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("Original vs Smoothed Expression: Top Variable Genes", 
                                           fontface = "bold", size = 14),
  p5c, ncol = 1, rel_heights = c(0.08, 0.92)
)

save_plot(p5c, file.path(plots_dir, "08a_original_vs_smoothed_top_genes"), 
                      width = 10, height = 7)
cat("  ✓ Saved original vs smoothed top genes scatter (PNG + PDF)\n")

# PLOT 5d: Density comparison for expression distributions
# Sample 3 genes: high, medium, low variance
sorted_var <- sort(var_original, decreasing = TRUE)
gene_high <- names(sorted_var)[1]
gene_med <- names(sorted_var)[floor(length(sorted_var)/2)]
gene_low <- names(sorted_var)[length(sorted_var) - 10]  # Avoid zero-variance genes

density_df <- rbind(
  data.frame(
    Expression = as.numeric(rna_original_train[gene_high, ]),
    Gene = gene_high,
    Type = "Original"
  ),
  data.frame(
    Expression = as.numeric(smoothed_train$rna_log1p[gene_high, ]),
    Gene = gene_high,
    Type = "Smoothed"
  ),
  data.frame(
    Expression = as.numeric(rna_original_train[gene_med, ]),
    Gene = gene_med,
    Type = "Original"
  ),
  data.frame(
    Expression = as.numeric(smoothed_train$rna_log1p[gene_med, ]),
    Gene = gene_med,
    Type = "Smoothed"
  ),
  data.frame(
    Expression = as.numeric(rna_original_train[gene_low, ]),
    Gene = gene_low,
    Type = "Original"
  ),
  data.frame(
    Expression = as.numeric(smoothed_train$rna_log1p[gene_low, ]),
    Gene = gene_low,
    Type = "Smoothed"
  )
)

density_df$Gene <- factor(density_df$Gene, levels = c(gene_high, gene_med, gene_low),
                          labels = c(paste0(gene_high, "\n(High var)"),
                                     paste0(gene_med, "\n(Med var)"),
                                     paste0(gene_low, "\n(Low var)")))

p5d <- ggplot(density_df, aes(x = Expression, fill = Type, color = Type)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  facet_wrap(~Gene, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Original" = "#D55E00", "Smoothed" = "#0072B2")) +
  scale_color_manual(values = c("Original" = "#D55E00", "Smoothed" = "#0072B2")) +
  labs(title = "Expression Density: Original vs Smoothed",
       subtitle = "Smoothing reduces sparsity and noise, especially for high-variance genes",
       x = "log₁p(CPM)", y = "Density") +
  theme_pub() +
  theme(legend.position = "bottom")

save_plot(p5d, file.path(plots_dir, "08b_expression_density_comparison"), 
                      width = 11, height = 5)
cat("  ✓ Saved expression density comparison (PNG + PDF)\n")

# PLOT 5e: Zero proportion before/after smoothing
zero_prop_original <- rowMeans(rna_original_train[sample_genes, ] == 0)
zero_prop_smoothed <- rowMeans(smoothed_train$rna_log1p[sample_genes, ] == 0)

zero_df <- data.frame(
  Original = zero_prop_original,
  Smoothed = zero_prop_smoothed
)

p5e <- ggplot(zero_df, aes(x = Original, y = Smoothed)) +
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

save_plot(p5e, file.path(plots_dir, "08c_dropout_reduction"), 
                      width = 7, height = 6)
cat("  ✓ Saved dropout reduction plot (PNG + PDF)\n")

# PLOT 5f: Mean expression preservation check
mean_original <- rowMeans(rna_original_train[sample_genes, ])
mean_smoothed <- rowMeans(smoothed_train$rna_log1p[sample_genes, ])

mean_df <- data.frame(
  Original = mean_original,
  Smoothed = mean_smoothed
)

cor_mean <- cor(mean_original, mean_smoothed)

p5f <- ggplot(mean_df, aes(x = Original, y = Smoothed)) +
  geom_point(alpha = 0.3, size = 1, color = "#E69F00") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#0072B2", linewidth = 0.8) +
  labs(title = "Mean Expression Preservation",
       subtitle = sprintf("Correlation: r = %.4f (smoothing preserves mean structure)", cor_mean),
       x = "Original Mean Expression", y = "Smoothed Mean Expression") +
  theme_pub()

save_plot(p5f, file.path(plots_dir, "08d_mean_expression_preservation"), 
                      width = 7, height = 6)
cat("  ✓ Saved mean expression preservation plot (PNG + PDF)\n")

# PLOT 5g: Heatmap of top genes before/after (side by side)
cat("Creating expression heatmaps (original vs smoothed)...\n")

# Select top 50 variable genes
top50_genes <- names(sort(var_original, decreasing = TRUE))[1:50]

# Sample 200 cells for visualization
set.seed(seed)
sample_cells <- sample(smoothed_train$cells, min(200, length(smoothed_train$cells)))

# Get expression matrices
orig_mat <- as.matrix(rna_original_train[top50_genes, sample_cells])
smooth_mat <- as.matrix(smoothed_train$rna_log1p[top50_genes, sample_cells])

# Scale for heatmap
orig_scaled <- t(scale(t(orig_mat)))
smooth_scaled <- t(scale(t(smooth_mat)))

# Create heatmap data frames
orig_hm_df <- reshape2::melt(orig_scaled)
colnames(orig_hm_df) <- c("Gene", "Cell", "Expression")
orig_hm_df$Type <- "Original"

smooth_hm_df <- reshape2::melt(smooth_scaled)
colnames(smooth_hm_df) <- c("Gene", "Cell", "Expression")
smooth_hm_df$Type <- "Smoothed"

hm_df <- rbind(orig_hm_df, smooth_hm_df)
hm_df$Type <- factor(hm_df$Type, levels = c("Original", "Smoothed"))

# Clip extreme values for visualization
hm_df$Expression <- pmax(pmin(hm_df$Expression, 3), -3)

p5g <- ggplot(hm_df, aes(x = Cell, y = Gene, fill = Expression)) +
  geom_tile() +
  facet_wrap(~Type, ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, limits = c(-3, 3),
                       name = "Scaled\nExpression") +
  labs(title = "Expression Heatmap: Original vs Smoothed",
       subtitle = sprintf("Top 50 variable genes × %d cells (z-scored)", length(sample_cells)),
       x = "Cells", y = "Genes") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 5))

save_plot(p5g, file.path(plots_dir, "08e_expression_heatmap_comparison"), 
                      width = 14, height = 10)
cat("  ✓ Saved expression heatmap comparison (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 6: RNA vs ATAC Correlation in Joint Space
# ----------------------------------------------------------------------------
cat("Creating RNA-ATAC correlation plots...\n")

# Calculate correlation between RNA PCA and ATAC LSI (first few components)
n_comp_plot <- min(5, pca_dims)

rna_pca_train <- smoothed_train$rna_pca[, 1:n_comp_plot]
atac_lsi_train <- smoothed_train$atac_lsi[, 1:n_comp_plot]

cor_matrix <- cor(rna_pca_train, atac_lsi_train)

cor_df <- expand.grid(
  RNA_PC = paste0("PC", 1:n_comp_plot),
  ATAC_LSI = paste0("LSI", 1:n_comp_plot)
)
cor_df$Correlation <- as.vector(cor_matrix)

p6 <- ggplot(cor_df, aes(x = RNA_PC, y = ATAC_LSI, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, limits = c(-1, 1),
                       name = "Correlation") +
  labs(title = "RNA PCA vs ATAC LSI Correlation",
       subtitle = "Cross-modality relationships in joint space",
       x = "RNA Principal Component", y = "ATAC LSI Component") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p6, file.path(plots_dir, "09_rna_atac_correlation_heatmap"), 
                      width = 7, height = 6)
cat("  ✓ Saved RNA-ATAC correlation heatmap (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 7: Summary Parameters Panel
# ----------------------------------------------------------------------------
cat("Creating summary panel...\n")

summary_text <- paste0(
  "METACELL CREATION PARAMETERS\n\n",
  "Dimensionality Reduction:\n",
  sprintf("  - RNA PCA dimensions: %d\n", pca_dims),
  sprintf("  - ATAC LSI dimensions: %d\n", pca_dims),
  sprintf("  - Variable features: %d\n\n", n_variable_features),
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

p7 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = summary_text, 
           hjust = 0.5, vjust = 0.5, size = 4, family = "mono") +
  theme_void() +
  theme(plot.margin = margin(20, 20, 20, 20),
        plot.background = element_rect(fill = "white", color = NA))

save_plot(p7, file.path(plots_dir, "10_summary_parameters"), 
                      width = 6, height = 5)
cat("  ✓ Saved summary parameters (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# COMBINED PDF REPORT (Vector Graphics)
# ----------------------------------------------------------------------------
cat("\nCreating combined PDF report...\n")

cairo_pdf(file.path(plots_dir, "Step_030_Metacell_Report.pdf"), width = 10, height = 8)

# Title page
plot.new()
text(0.5, 0.6, "Step 030: Metacell Creation Report", cex = 2, font = 2)
text(0.5, 0.45, sprintf("Sample: %s", SAMPLE_NAME), cex = 1.4)
text(0.5, 0.35, sprintf("Date: %s", Sys.Date()), cex = 1.2)
text(0.5, 0.25, "Ready Figures", cex = 1, font = 3)

# Print plots
if (exists("p1a")) print(p1a)
if (exists("p1b")) print(p1b)
if (exists("p2a")) print(p2a)
if (exists("p2b")) print(p2b)
print(p3)
print(p4)
print(p5a)
print(p5b)
if (exists("p5c")) print(p5c)
if (exists("p5d")) print(p5d)
if (exists("p5e")) print(p5e)
if (exists("p5f")) print(p5f)
if (exists("p5g")) print(p5g)
print(p6)
print(p7)

dev.off()
cat("  ✓ Saved combined PDF report (vector graphics)\n")

cat("\n")
cat("plots saved to:", plots_dir, "\n")
cat("  - Each plot saved as both .png (600 dpi) and .pdf (vector)\n")
cat("  - Combined report: Step_030_Metacell_Report.pdf\n")

# ============================================================================
# COMPLETION MESSAGE (FINAL)
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 030 COMPLETE\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("k-NN smoothing completed for all data splits\n\n")

cat("Parameters:\n")
cat(sprintf("  - k = %d neighbors\n", k_neighbors))
cat(sprintf("  - PCA/LSI dims = %d\n", pca_dims))
cat(sprintf("  - Variable genes = %d\n", n_variable_features))
cat(sprintf("  - Random seed = %d\n\n", seed))

cat("Key features:\n")
cat("  ✓ Training: Built PCA/LSI transformations from scratch\n")
cat("  ✓ Val/Test: Projected into training coordinate system\n")
cat("  ✓ Val/Test: Smoothed only among val/test cells (no data leakage)\n")
cat("  ✓ All cells preserved (1-to-1 mapping, no aggregation)\n")
cat("  ✓ plots (PNG + PDF)\n\n")

cat("Output files:\n")
cat("  - train_transforms.rds (PCA/LSI/z-score parameters)\n")
cat("  - smoothed_train.rds (smoothed training data + coordinates)\n")
cat("  - smoothed_validation.rds (smoothed validation data + coordinates)\n")
cat("  - smoothed_test.rds (smoothed test data + coordinates)\n\n")

cat("Diagnostic plots (in plots/ subdirectory, PNG + PDF formats):\n")
cat("  01 - UMAP of joint space by split\n")
cat("  02 - UMAP faceted by split\n")
cat("  03 - PCA variance explained\n")
cat("  04 - LSI variance explained\n")
cat("  05 - Split distribution bar chart\n")
cat("  06 - k-NN distance distribution\n")
cat("  07 - Smoothing variance comparison (original vs smoothed)\n")
cat("  08 - Variance reduction histogram\n")
cat("  08a - Original vs smoothed: top variable genes scatter\n")
cat("  08b - Expression density comparison (original vs smoothed)\n")
cat("  08c - Dropout reduction plot\n")
cat("  08d - Mean expression preservation\n")
cat("  08e - Expression heatmap comparison\n")
cat("  09 - RNA-ATAC correlation heatmap\n")
cat("  10 - Summary parameters\n")
cat("  + Combined PDF report (vector graphics)\n\n")

cat("Output directory:", OUTPUT_METACELLS_DIR, "\n\n")

cat("Next step: Run 04_feature_extraction.R\n")
