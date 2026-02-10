# Step 030: Automated Metacell Creation via KNN Smoothing (ATAC-ONLY)
# Automated Metacell Creation via KNN Smoothing (ATAC-ONLY, LSI-based)
# Purpose: Denoise single-cell ATAC data by smoothing in LSI space using
#          cell-wise k-NN weights computed within each split (train/val/test).
#
# CRITICAL DESIGN NOTES:
#  - Training LSI is fit on TRAIN only; validation/test are projected into
#    the training coordinate system using stored IDF and loadings (prevents leakage).
#  - The first LSI component is dropped (captures sequencing depth).
#  - Smoothing is performed within-split (no mixing between train/val/test).
#  - This automated script reads configuration from `config.R` and expects
#    `INPUT_SEURAT_SPLITS` to point to a Seurat RDS that already contains
#    the `data_split` assignment.
#
# OUTPUTS (saved to `ATAC_METACELL_OUTPUT_DIR`)
#  - train_transforms.rds : list(L, idf, training_peaks, z_mu, z_sd, lsi_dims, k_neighbors)
#  - smoothed_train.rds / smoothed_validation.rds / smoothed_test.rds : lists containing
#    smoothed ATAC (log1p CPM), raw smoothed counts, LSI embeddings, RNA targets (if present).
#
# USAGE:
#  - Edit `config.R` to set `INPUT_SEURAT_SPLITS` and `ATAC_METACELL_OUTPUT_DIR`.
#  - Run: `Rscript 03a_metacell_creation.R`
#
# References: Signac/Seurat functions used: FindTopFeatures, RunTFIDF, RunSVD/RunLSI

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Signac)
  library(Matrix)
  library(RANN)
  library(dplyr)
})

# ---------------- Load Config ----------------
if (!exists("CONFIG_LOADED")) {
  config_file <- "config.R"
  if (!file.exists(config_file)) stop("ERROR: config.R not found in working directory.")
  source(config_file)
}

# ---------------- Load Config ----------------
# Load pipeline-wide configuration. The `config.R` should set paths like
# `INPUT_SEURAT_SPLITS` (multiome-derived Seurat with splits) and output
# directories. The script accepts both `ATAC_*` aliases and legacy
# `OUTPUT_*` names; these aliases are resolved below.
if (!exists("ATAC_INPUT_SPLIT_DIR")) {
  if (exists("OUTPUT_SPLITS_DIR")) {
    ATAC_INPUT_SPLIT_DIR <- dirname(path.expand(OUTPUT_SPLITS_DIR))
  } else {
    ATAC_INPUT_SPLIT_DIR <- path.expand(".")
  }
}
if (!exists("ATAC_METACELL_OUTPUT_DIR")) {
  if (exists("OUTPUT_METACELLS_DIR")) {
    ATAC_METACELL_OUTPUT_DIR <- path.expand(OUTPUT_METACELLS_DIR)
  } else if (exists("OUTPUT_METACELLS_DIR")) {
    ATAC_METACELL_OUTPUT_DIR <- path.expand(OUTPUT_METACELLS_DIR)
  } else {
    ATAC_METACELL_OUTPUT_DIR <- file.path(path.expand("."), "metacells")
  }
}
if (!exists("ATAC_K_NEIGHBORS") && exists("K_NEIGHBORS")) ATAC_K_NEIGHBORS <- K_NEIGHBORS
if (!exists("ATAC_LSI_DIMS") && exists("LSI_DIMS")) ATAC_LSI_DIMS <- LSI_DIMS
if (!exists("ATAC_RANDOM_SEED") && exists("SEED_METACELL")) ATAC_RANDOM_SEED <- SEED_METACELL
if (!exists("ATAC_MIN_PEAK_QUANTILE")) ATAC_MIN_PEAK_QUANTILE <- "q5"

## Backwards-compatible alias: some older scripts expect `output_dir`
if (!exists("output_dir")) output_dir <- ATAC_METACELL_OUTPUT_DIR
output_dir <- path.expand(output_dir)
ATAC_METACELL_OUTPUT_DIR <- path.expand(ATAC_METACELL_OUTPUT_DIR)

cat("Input directory:", ATAC_INPUT_SPLIT_DIR, "\n")
cat("Output directory:", output_dir, "\n")
if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE); cat("Created output directory\n") }

k_neighbors <- ATAC_K_NEIGHBORS
lsi_dims    <- ATAC_LSI_DIMS
seed        <- ATAC_RANDOM_SEED
min_cutoff  <- ATAC_MIN_PEAK_QUANTILE
set.seed(seed)

# ---------------- Load Data ----------------
## Purpose: load the Seurat object that already contains train/validation/test
## split assignments. This script requires an explicit multiome-derived
## Seurat splits RDS path set in `config.R` as `INPUT_SEURAT_SPLITS`.


if (!exists("INPUT_SEURAT_SPLITS") || !nzchar(INPUT_SEURAT_SPLITS)) {
  stop("CONFIGURATION ERROR: `INPUT_SEURAT_SPLITS` must be set in config.R to the Seurat splits RDS produced by the multiome pipeline.")
}
input_file <- path.expand(INPUT_SEURAT_SPLITS)
if (!file.exists(input_file)) {
  stop(sprintf("Seurat splits RDS not found at INPUT_SEURAT_SPLITS: %s\nPlease run the multiome pipeline or set INPUT_SEURAT_SPLITS correctly in config.R.", input_file))
}
cat("Loading Seurat splits from:", input_file, "\n")

seurat_obj <- readRDS(input_file)
cat(sprintf("Loaded object: %d cells, %d peaks\n", ncol(seurat_obj), nrow(seurat_obj[["ATAC"]])))

# ---------------- Input validation ----------------
# Ensure the Seurat object contains the expected assays and split assignment.
# This guards against silent failures later in the pipeline and provides
# clear, actionable errors for the user.
if (!"ATAC" %in% names(seurat_obj@assays)) {
  stop("ERROR: ATAC assay not found in Seurat object")
}
if (!"data_split" %in% colnames(seurat_obj@meta.data)) {
  stop("ERROR: data_split column not found in Seurat object metadata. Ensure INPUT_SEURAT_SPLITS contains a 'data_split' column with values 'train','validation','test'.")
}

# Report per-split cell counts and simple sanity checks
for (split in c("train", "validation", "test")) {
  n_cells <- sum(seurat_obj$data_split == split)
  cat(sprintf("  %s: %d cells\n", split, n_cells))
  if (n_cells == 0) stop(sprintf("ERROR: No cells found for %s split. Check data_split assignments.", split))
  if (n_cells <= k_neighbors) warning(sprintf("WARNING: %s has %d cells (≤ k=%d). k will be adjusted.", split, n_cells, k_neighbors))
}
cat("Input validation passed.\n")

# ---------------- Helper Functions ----------------
# Zscore calculation
z <- function(M) {
  out <- scale(M, center = TRUE, scale = TRUE)
  if (is.matrix(out)) rownames(out) <- rownames(M)
  out
}

#-----------------------------------------------------------------------------------------
# build_knn_smooth
# Purpose:
#     Build a cell x cell smoothing matrix A for k-NN smoothing
#     When you multiply:
#       smoothed_counts = atac_counts %*% A
#     basically each cell's counts become the average of its own counts + its neighbors(k)
#     nearest neighbors in the given space
#
# Arugments:
#     scores : numeric matrix , rows = cells, cols = dimensions (LSI)
#     cell_names : character vector of length nrow(scores) : cell barcodes
#     k_smooth = number of neighbors (including self) to use per cell
#     seed = random seed (for reproducibility)

# Returns:
#     A : sparse matrix n x n matrix (dgCMatrix), with:
#         - rows = source cells(neighbors j)
#         - columns = target cells (centers i)
#       Each column sums to 
#------------------------------------------------------------------------------------------
build_knn_smooth <- function(scores, cell_names, k_smooth = 20, seed = 1) {
  set.seed(seed)
  n <- nrow(scores)
  # cap k so it never exceeds number of cells - 1 (neighbors excluding self)
  total_k <- k_smooth
  k_use <- max(1, min(total_k - 1, n - 1))
  if (n <= 1) {
    warning(sprintf("Too few cells (%d) for k-NN smoothing. Using identity.", n))
    return(sparseMatrix(
      i = 1, j = 1, x = 1,
      dims = c(n, n),
      dimnames = list(cell_names, cell_names)
    ))
  }
  # Get neighbors (self included) 
  # For each cell, find the nearest neighbors based on Euclidean distance 
  nn <- RANN::nn2(
    data       = scores,
    query      = scores,
    k          = k_use + 1,     # self + k neighbors
    searchtype = "standard"
  )
  # idx = n × (k_use+1), including self
  idx <- nn$nn.idx
  # Build smoothing matrix
  total_k_eff <- k_use + 1
  weight  <- 1 / total_k_eff
  
  # Build indices: for each center cell i, for each neighbor j in idx[i, ]
  i_idx <- as.vector(t(idx))                 # row indices (neighbor j)
  j_idx <- rep(seq_len(n), each = total_k)   # col indices (center cell i)
  x_vals <- rep(weight, length(i_idx))
  
  A <- sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = x_vals,
    dims = c(n, n),
    dimnames = list(cell_names, cell_names)
  )
  
  stopifnot(
    identical(rownames(A), cell_names),
    identical(colnames(A), cell_names)
  )
  A
}

# --------------------------------------------------------------------------------------------
# Apple k-NN smoothing to raw counts and normalize
#     Each column of A_smooth contains the weights for smoothing a single cell
#     So: smoothed_counts = raw_counts %*% A_smooth
#     Matrix multiplication: peaks/genes x cells
# -------------------------------------------------------------------------------------------
smooth_atac_accessibility <- function(atac_counts, A_smooth) {
  
  # enforce exact order: columns of counts == columns/rows of A
  if (!identical(colnames(atac_counts), colnames(A_smooth))) {
    atac_counts <- atac_counts[, colnames(A_smooth), drop = FALSE]
  }
  if (!identical(rownames(A_smooth), colnames(A_smooth))) {
    A_smooth <- A_smooth[colnames(A_smooth), colnames(A_smooth)]
  }
  
  # Library size Normalization for each smoothed cell
  # a. Library size : total counts per cell (after smoothing)
  # b. Convert smoothed counts to CPM
  # CPM = counts per million
  # c. Log-transform (log1p)
  
  atac_smoothed <- atac_counts %*% A_smooth
  atac_lib <- Matrix::colSums(atac_smoothed)
  atac_cpm <- t(t(atac_smoothed) / (atac_lib + 1e-8)) * 1e6
  atac_log1p <- log1p(atac_cpm)                        
  colnames(atac_log1p) <- colnames(A_smooth)
  list(atac_log1p = atac_log1p, atac_lib = atac_lib, atac_counts = atac_smoothed)
}

# -------------------------------------------------------------------------------------------
# Build TRAIN transforms (fit) 
#   - Subset to training cells
#   - Run TF-IDF + LSI on Train Only
#   - Drop depth component
#   - Z-score LSI dimensions using TRAIN statistics
#   - Build k-NN smootthin matrix 
#   - Apply smoothing to raw ATAC counts
#-------------------------------------------------------------------------------------------

cat("\n=== Processing TRAINING data (ATAC-ONLY) ===\n")

# 1. Subset Seurat obj to TRAIN cells only
train_cells  <- colnames(seurat_obj)[seurat_obj$data_split == "train"]
seurat_train <- seurat_obj[, train_cells]

# Set Assay to ATAC
DefaultAssay(seurat_train) <- "ATAC"

# 2. Ensure ATAC counts are stored as a sparse counts
ac <- GetAssayData(seurat_train, assay = "ATAC", layer = "counts")
if (!inherits(ac, "dgCMatrix")) {
  seurat_train <- SetAssayData(
    seurat_train, 
    assay = "ATAC", 
    layer = "counts", 
    new.data = as(ac, "dgCMatrix"))
}
cat(sprintf("Peaks: %d\n", nrow(seurat_train[["ATAC"]])))

# 3. Select top ATAC peaks and compute TF-IDF 
#   - min.cutoff <- "q5" keeps peaks above 5th percentile 
seurat_train <- FindTopFeatures(seurat_train, min.cutoff = "q5")
seurat_train <- RunTFIDF(seurat_train, verbose = FALSE)

# 4.Run LSI (SVD) on TRAIN TF-IDF to get low-dimensional representation
seurat_train <- RunSVD(seurat_train, n = lsi_dims, verbose = FALSE)

# Extract TRAIN LSI embeddings (cells x components) and drop first component
# LSI1 usually captured sequencing depth
lsi_train <- Embeddings(seurat_train, "lsi")[, -1, drop = FALSE]
rownames(lsi_train) <- train_cells

cat(sprintf("ATAC LSI (train): %d cells × %d comps (dropped 1st)\n",
            nrow(lsi_train), ncol(lsi_train)))

# 5. Define helper function for z-scoring columns of a matrix
# z_fit_cols: compute per-dimension and sd
z_fit_cols <- function(M){
  mu <- colMeans(M)
  sd <- apply(M, 2, sd)
  # Avoid division by zero
  sd[!is.finite(sd) | sd < 1e-8] <- 1
  list(mu = mu, sd = sd)
}

#   z_apply: apply z_scoring given mu and sd
z_apply_cols <- function(M, mu, sd){
  sweep(sweep(M, 2, mu, "-"), 2, sd, "/")
}

# Fit z_score parameters on TRAIN LSI only
z_params <- z_fit_cols(lsi_train)

# Apply z_scoring to TRAIN LSI (JE_train = joint embedding (standardized))
JE_train <- z_apply_cols(lsi_train, z_params$mu, z_params$sd)

# 6. Build k-NN smoothing matrix on TRAIN cells in z-scored LSI space
#     - Each column i of A_smooth_train contains weights of neighbors j
#       contributing to smoothed cell i
A_smooth_train <- build_knn_smooth(
  scores = JE_train, 
  cell_names = train_cells, 
  k_smooth = k_neighbors,
  seed = seed)

# 7. Smooth raw ATAC counts for TRAIN cells using A_smooth_train
atac_counts_train <- GetAssayData(seurat_train, assay = "ATAC", layer = "counts")

# smooth_atac_accessibility:
#   - ensures peak/cell order matches A_smooth_train
#   - computes counts %*% A
#   - normalizes and log1p transforms
smoothed_train <- smooth_atac_accessibility(atac_counts_train, A_smooth_train)


# ---------------------------------------------------------------------------------------------
# Compute and save TRAIN LSI transformations (IDF + loadings)
# These will be used to project to the validation/test ATAC into
# training's LSI space without refitting
#----------------------------------------------------------------------------------------------

# Set assay to ATAC
DefaultAssay(seurat_train) <- "ATAC"

# 1. Extract raw peak x cell count matrix
counts_train <- GetAssayData(seurat_train, assay = "ATAC", layer = "counts")
if (!inherits(counts_train, "dgCMatrix")) counts_train <- as(counts_train, "dgCMatrix")

# counts_train dimensions:
#   rows = peaks
#   columns = training cells

# 2. Extract LSI Loadings from learned on training set
#     Loadings : peaks x LSI dimensions (before dropping component 1)
L <- Loadings(seurat_train[["lsi"]])             # peaks × dims (includes comp1)
training_peaks <- rownames(L)       # store exact peak order

# Drop comp1 in loadings (sequencing depth)
L <- L[, -1, drop = FALSE]


# 3. Compute IDF weights for TF-IDF:
#     IDF = log(1+(#cells / #cells-with-peak))
n_cells_train <- ncol(counts_train)                     # number of training cells
nz_per_peak   <- Matrix::rowSums(counts_train > 0)      # #cells where peak is accessible
idf <- log1p(n_cells_train / pmax(1, nz_per_peak))
names(idf) <- rownames(counts_train)                    # name vector for alignment

train_transforms <- list(
  L              = L,                # loadings (peaks × dims, comp1 removed)
  idf            = idf,              # training IDF per peak
  training_peaks = training_peaks,   # order anchor
  z_mu           = z_params$mu,      # z-score params
  z_sd           = z_params$sd,
  lsi_dims       = lsi_dims,
  k_neighbors    = k_neighbors
)
saveRDS(train_transforms, file.path(output_dir, "train_transforms.rds"))
cat("[OK] Saved train_transforms.rds\n")

# Add this before saving train_output:
DefaultAssay(seurat_train) <- "RNA"
seurat_train <- NormalizeData(seurat_train, verbose = FALSE)
rna_data_train <- GetAssayData(seurat_train, layer = "data")

# Save smoothed TRAIN bundle
train_output <- list(
  atac_log1p = smoothed_train$atac_log1p,
  atac_lib   = smoothed_train$atac_lib,
  atac_counts= smoothed_train$atac_counts,
  atac_lsi   = lsi_train,
  cells      = train_cells,
  rna_log1p  = rna_data_train, # unsmoothed rna
  split_name = "train",
  processing_seed = seed,
  k_smooth   = k_neighbors,
  lsi_dims   = lsi_dims
)
saveRDS(train_output, file.path(output_dir, "smoothed_train.rds"))
cat("[OK] Saved smoothed_train.rds (with RNA)\n")

# ---------------------------------------------------------------------------------------------
# Project VAL/TEST split into TRAIN LSI space
# using TRAIN Loadings and TRAIN IDF (no refit)
#
# Arguments:
#     seurat_split: Seurat object containing only val or test cells
#     L           :TRAIN LSI Loadings (peaks x dims, component 1 dropped)
#     idf         : TRAIN IDF vector for peaks
#     training_peaks : charactor vector of TRAIN peak_names (order anchor)
#     assay: which assay to use
#
# Returns:
#     emb: matrix (#cells_in_split x #dims), LSI embeddings in TRAIN space
#-----------------------------------------------------------------------------------------------
cat("\n=== Processing VALIDATION and TEST data (ATAC-ONLY) ===\n")

# Load Train transformations
train_transforms <- readRDS(file.path(output_dir, "train_transforms.rds"))

# Helper function for projection
project_split_lsi <- function(seurat_split, L, idf, training_peaks, assay = "ATAC") {
  # Set assay to ATAC
  DefaultAssay(seurat_split) <- assay

  # 1. Extract raw ATAC counts (peaks × cells) for this split
  Xc <- GetAssayData(seurat_split, assay = assay, layer = "counts")
  if (!inherits(Xc, "dgCMatrix")) Xc <- as(Xc, "dgCMatrix")

  # 2. Align peaks to TRAIN peaks:
  #     - Find intersection between TRAIN peaks and split peaks
  #     - Keep TRAIN order to ensure consistent alignment
  common <- training_peaks[training_peaks %in% rownames(Xc)]
  
  # Safety check: ensure enough overlapping peaks to make projection work well
  if (length(common) < max(500, floor(0.5 * length(training_peaks)))) {
    stop(sprintf("Too few overlapping peaks for projection: %d/%d",
                 length(common), length(training_peaks)))
  }
  
  # Subset and reorder the split counts, loadings and IDF to common peaks
  Xc   <- Xc[common, , drop = FALSE]    # peaks x cells
  Luse <- L[common, , drop = FALSE]     # peaks x dims
  idf  <- idf[common]                   # IDF vector for these peaks

  # 3. Compute TF(term frequency) (per cell) for this split:
  #   TF[p, i] = counts[p, i] / total_counts[i]
  lib <- Matrix::colSums(Xc) 
  lib[lib == 0] <- 1                    # avoid division by zero
  
  # t(Xc) / lib : cells x peaks, then transpose back to peaks x cells
  TF  <- t(t(Xc) / lib)                 # peaks × cells

  # 4. Apply training IDF to get TF-IDF
  #       - No recomputing IDF on val/test (prevent data leakage and ensures
  #         the same scaling of rare/common peaks as in training)
  TFIDF <- TF * idf     # broadcasts idf over columns (peaks x cells)
  
  # 5. Project into TRAIN LSI space:
  #    emb = (cells x peaks) %*% (peaks x dims) = cells x dims
  emb   <- as.matrix(Matrix::t(TFIDF) %*% Luse)

  # Label rows by cell barcodes
  rownames(emb) <- colnames(seurat_split)
  
  # Return LSI embeddings of this split in TRAIN coordinate system
  emb
}


for (split_name in c("validation", "test")) {
  cat(sprintf("\n--- Processing %s ---\n", toupper(split_name)))

  split_cells  <- colnames(seurat_obj)[seurat_obj$data_split == split_name]
  seurat_split <- seurat_obj[, split_cells]

  # ensure sparse
  DefaultAssay(seurat_split) <- "ATAC"
  ac <- GetAssayData(seurat_split, assay = "ATAC", slot = "counts")
  if (!inherits(ac, "dgCMatrix")) {
    seurat_split <- SetAssayData(seurat_split, assay = "ATAC", slot = "counts",
                                 new.data = as(ac, "dgCMatrix"))
  }

  # Project into training LSI (no refit)
  lsi_split <- project_split_lsi(
    seurat_split,
    L              = train_transforms$L,
    idf            = train_transforms$idf,
    training_peaks = train_transforms$training_peaks,
    assay = "ATAC"
  )

  # Z-score with TRAIN μ/σ (no leakage)
  JE_split <- z_apply_cols(lsi_split, train_transforms$z_mu, train_transforms$z_sd)

  # kNN within split only
  A_smooth_split <- build_knn_smooth(JE_split, split_cells, k_smooth = k_neighbors, seed = seed)

  # Smooth counts
  atac_counts_split <- GetAssayData(seurat_split, assay = "ATAC", slot = "counts")
  smoothed_split <- smooth_atac_accessibility(atac_counts_split, A_smooth_split)
  cat(sprintf("[OK] %s: %d cells smoothed\n", split_name, ncol(smoothed_split$atac_log1p)))

  # RNA targets for metrics
  DefaultAssay(seurat_split) <- "RNA"
  seurat_split <- NormalizeData(seurat_split, verbose = FALSE)
  rna_data_split <- GetAssayData(seurat_split, slot = "data")

  # Save
  split_output <- list(
    atac_log1p = smoothed_split$atac_log1p,
    atac_lib   = smoothed_split$atac_lib,
    atac_counts= smoothed_split$atac_counts,
    atac_lsi   = lsi_split,
    cells      = split_cells,
    rna_log1p  = rna_data_split,
    split_name = split_name,
    processing_seed = seed,
    k_smooth   = k_neighbors,
    lsi_dims   = lsi_dims
  )
  saveRDS(split_output, file.path(output_dir, paste0("smoothed_", split_name, ".rds")))
  cat(sprintf("[OK] Saved smoothed_%s.rds (with RNA)\n", split_name))
}


cat("\n==========================\n")
cat("Step 030 Complete (ATAC-ONLY)\n")
cat("==========================\n")
cat(sprintf("k=%d | LSI dims=%d | seed=%d\n", k_neighbors, lsi_dims, seed))


# ===============Summary=========================
cat("\n" , paste(rep("=", 60), collapse = ""), "\n")
cat("=== Step 030 Complete (ATAC-ONLY) ===\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")
cat("k-NN smoothing completed for all data splits (ATAC-ONLY)\n")
cat(sprintf("Parameters:\n"))
cat(sprintf("  - k = %d neighbors\n", k_neighbors))
cat(sprintf("  - LSI dims = %d\n", lsi_dims))
cat(sprintf("  - Random seed = %d\n", seed))
cat("\nKey features:\n")
cat("  ✓ Training: Built LSI transformation from scratch using Signac\n")
cat("  ✓ Val/Test: Projected into training coordinate system\n")
cat("  ✓ Val/Test: Smoothed only among val/test cells (no data leakage)\n")
cat("  ✓ All cells preserved (1-to-1 mapping, no aggregation)\n")
cat("  ✓ ATAC-ONLY modeling: Only ATAC peaks used as features\n")
cat("  ✓ RNA expression included as target variable (y) for supervised learning\n")
cat("  ✓ LSI: Skipped 1st component (sequencing depth)\n")
cat("  ✓ IDF values properly stored and reused for projection\n")
cat("  ✓ Input validation: Checked assays, splits, and cell counts\n")
cat("\nFiles saved:\n")
cat("  - train_transforms.rds (LSI parameters with IDF)\n")
cat("  - smoothed_train.rds (smoothed ATAC + RNA expression)\n")
cat("  - smoothed_validation.rds (smoothed ATAC + RNA expression)\n")
cat("  - smoothed_test.rds (smoothed ATAC + RNA expression)\n")
cat("\nNext: Use smoothed ATAC as features (X) and RNA as targets (y) for modeling\n")

