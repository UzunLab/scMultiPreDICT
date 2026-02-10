# ============================================================================
# Step_030: Metacell Creation via k-NN Smoothing (scRNA-only, Automated)
# ============================================================================
# This script creates metacells using k-NN smoothing in RNA PCA space.
# Uses RNA data only from the multiome Seurat object (ignores ATAC).
#
# Input: Seurat object with data splits (from multiome pipeline)
# Output: Smoothed RNA expression matrices for train/val/test splits
#
# Key Features:
#   - PCA built from TRAINING data only
#   - Val/Test projected into training PCA space
#   - Smoothing done WITHIN each split (no data leakage)
#   - All cells preserved (1-to-1 mapping)
#   - Publication-ready diagnostic plots
#
# Usage:
#   1. Edit config.R with your parameters
#   2. Run: Rscript 03a_metacell_creation.R
# ============================================================================

if (!exists("CONFIG_LOADED")) {
  # Get the directory of this script to source config.R
  get_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      return(dirname(normalizePath(sub("^--file=", "", file_arg))))
    }
    return(".")
  }
  script_dir <- get_script_dir()

  # Source configuration file
  config_path <- file.path(script_dir, "config.R")
  if (!file.exists(config_path)) {
    config_path <- "config.R"
  }

  if (!file.exists(config_path)) {
    stop("config.R not found! Please ensure config.R is in the same directory as this script.")
  }

  cat("Loading configuration from:", config_path, "\n")
  source(config_path)
}

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(RANN)         # Fast k-NN search
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(Cairo)
  library(viridis)
})

# Ensure irlba is available
if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("Package 'irlba' is required. Install with: install.packages('irlba')")
}

# ============================================================================
# PUBLICATION-READY PLOTTING SETUP
# ============================================================================

# Colorblind-friendly palette (Wong 2011)
colorblind_palette <- c(
  "#E69F00",  # Orange
  "#56B4E9",  # Sky blue
  "#009E73",  # Bluish green
  "#F0E442",  # Yellow
  "#0072B2",  # Blue
  "#D55E00",  # Vermillion
  "#CC79A7",  # Reddish purple
  "#999999"   # Gray
)

# Split colors for publication
split_colors_pub <- c(
  "train" = "#0072B2",       # Blue
  "validation" = "#E69F00",  # Orange
  "test" = "#009E73"         # Green
)

# Plotting theme
theme_pub <- function(base_size = 12, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5, 
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 2, hjust = 1, color = "gray50"),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
      strip.text = element_text(size = base_size, face = "bold", color = "black"),
      plot.margin = margin(15, 15, 15, 15)
    )
}

# Save plot function
save_plot <- function(plot, filename, width = 10, height = 8, dpi = 600, 
                                   create_pdf = TRUE) {
  ggsave(
    filename = paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  cat(sprintf("  Saved: %s.png\n", basename(filename)))
  
  if (create_pdf) {
    ggsave(
      filename = paste0(filename, ".pdf"),
      plot = plot,
      width = width,
      height = height,
      device = cairo_pdf,
      bg = "white"
    )
    cat(sprintf("  Saved: %s.pdf\n", basename(filename)))
  }
}

# ============================================================================
# PRINT CONFIGURATION
# ============================================================================
print_config()
print_output_directories()

# Create output directories
create_output_directories()

# Create plots directory
plots_dir <- file.path(OUTPUT_METACELLS_DIR, "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 050: Metacell Creation via k-NN Smoothing (scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat(sprintf("Dimensionality Reduction Method: %s\n", DIM_REDUCTION_METHOD))
cat(sprintf("Output directory: %s\n\n", OUTPUT_METACELLS_DIR))

# ============================================================================
# RNA PCA FUNCTIONS
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
    train_scores = pr$x,
    sdev = pr$sdev
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

# ============================================================================
# K-NN SMOOTHING FUNCTIONS
# ============================================================================

# Build cell x cell smoothing matrix A for k-NN smoothing
build_knn_smooth <- function(scores, cell_names, k_smooth = 20, seed = 1) {
  set.seed(seed)
  n <- nrow(scores)
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
  
  nn <- RANN::nn2(
    data       = scores,
    query      = scores,
    k          = k_use + 1,
    searchtype = "standard"
  )
  
  idx <- nn$nn.idx
  total_k_eff <- k_use + 1
  weight  <- 1 / total_k_eff
  
  i_idx <- as.vector(t(idx))
  j_idx <- rep(seq_len(n), each = total_k)
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

# Smooth RNA expression
smooth_rna_expression <- function(rna_counts, A_smooth) {
  rna_smoothed <- rna_counts %*% A_smooth
  
  rna_lib <- Matrix::colSums(rna_smoothed)
  rna_cpm <- t(t(rna_smoothed) / (rna_lib + 1e-8)) * 1e6
  rna_log1p <- log1p(rna_cpm)
  
  list(
    rna_log1p = rna_log1p,
    rna_lib = rna_lib,
    rna_counts = rna_smoothed
  )
}

# ============================================================================
# LOAD INPUT DATA
# ============================================================================
cat("=== Loading Seurat object with data splits ===\n\n")

input_file <- path.expand(INPUT_SEURAT_SPLITS)
cat("Input file:", input_file, "\n")

if (!file.exists(input_file)) {
  stop("ERROR: Input file not found: ", input_file, "\n",
       "Please ensure the multiome preprocessing pipeline was run first.")
}

seurat_obj <- readRDS(input_file)

cat(sprintf("Loaded object: %d cells, %d genes (RNA)\n", 
            ncol(seurat_obj), nrow(seurat_obj[["RNA"]])))

# Check for data splits
if (!"data_split" %in% colnames(seurat_obj@meta.data)) {
  stop("ERROR: 'data_split' column not found in metadata. Run Step_040 first.")
}

split_counts <- table(seurat_obj$data_split)
cat("Data splits:\n")
for (split in names(split_counts)) {
  cat(sprintf("  %s: %d cells\n", split, split_counts[split]))
}

# Set random seed
set.seed(SEED_METACELL)

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
  seurat_train <- SetAssayData(seurat_train, assay = "RNA", layer = "counts", 
                               new.data = as(rc, "dgCMatrix"))
}

cat(sprintf("Genes: %d\n", nrow(seurat_train[["RNA"]])))

# Build RNA PCA from training data
cat(sprintf("\n[INFO] Building RNA PCA (dims=%d, features=%d)...\n", 
            PCA_DIMS, N_VARIABLE_FEATURES))
rna_fit <- build_rna_pca_fit(seurat_train, pca_dims = PCA_DIMS, 
                              nfeatures = N_VARIABLE_FEATURES)
rna_scores_train <- rna_fit$train_scores
rownames(rna_scores_train) <- train_cells

cat(sprintf("  Variable genes: %d\n", length(rna_fit$vf)))
cat(sprintf("  PCA components: %d\n", ncol(rna_scores_train)))

# Build k-NN smoothing matrix for training
cat(sprintf("\n[INFO] Building k-NN smoothing matrix (k=%d)...\n", K_NEIGHBORS))
A_smooth_train <- build_knn_smooth(rna_scores_train, train_cells, 
                                    k_smooth = K_NEIGHBORS, seed = SEED_METACELL)

neighbors_per_cell <- Matrix::rowSums(A_smooth_train > 0)
cat(sprintf("  k-NN stats: mean=%.1f, min=%d, max=%d neighbors\n",
            mean(neighbors_per_cell), min(neighbors_per_cell), max(neighbors_per_cell)))

# Smooth training expression
cat("[INFO] Smoothing training expression...\n")
rna_counts_train <- GetAssayData(seurat_train, assay = "RNA", layer = "counts")
smoothed_train <- smooth_rna_expression(rna_counts_train, A_smooth_train)
cat(sprintf("[OK] Training: %d cells smoothed\n", ncol(smoothed_train$rna_log1p)))

# Save training transformation parameters
train_transforms <- list(
  rna_fit = rna_fit,
  pca_dims = PCA_DIMS,
  k_neighbors = K_NEIGHBORS,
  n_variable_features = N_VARIABLE_FEATURES
)
saveRDS(train_transforms, file.path(OUTPUT_METACELLS_DIR, "train_transforms.rds"))
cat("[OK] Saved train_transforms.rds\n")

# Save smoothed training data
train_output <- list(
  rna_log1p = smoothed_train$rna_log1p,
  rna_lib = smoothed_train$rna_lib,
  rna_counts = smoothed_train$rna_counts,
  rna_pca = rna_scores_train,
  cells = train_cells,
  split_name = "train",
  processing_seed = SEED_METACELL,
  k_smooth = K_NEIGHBORS,
  pca_dims = PCA_DIMS
)
saveRDS(train_output, file.path(OUTPUT_METACELLS_DIR, "smoothed_train.rds"))
cat("[OK] Saved smoothed_train.rds\n")

# Store for combined plots
all_pca_scores <- list(train = rna_scores_train)
all_splits_info <- list(train = list(n_cells = length(train_cells)))

# ============================================================================
# PROCESS VALIDATION AND TEST DATA
# ============================================================================
cat("\n=== Processing VALIDATION and TEST data ===\n\n")

for (split_name in c("validation", "test")) {
  cat(sprintf("\n--- Processing %s ---\n", toupper(split_name)))
  
  split_cells <- colnames(seurat_obj)[seurat_obj$data_split == split_name]
  seurat_split <- seurat_obj[, split_cells]
  
  cat(sprintf("%s cells: %d\n", split_name, length(split_cells)))
  
  # Normalize
  DefaultAssay(seurat_split) <- "RNA"
  seurat_split <- normalize_rna(seurat_split)
  
  rc <- GetAssayData(seurat_split, assay = "RNA", layer = "counts")
  if (!inherits(rc, "dgCMatrix")) {
    seurat_split <- SetAssayData(seurat_split, assay = "RNA", layer = "counts",
                                 new.data = as(rc, "dgCMatrix"))
  }
  
  # Project into training PCA space
  cat(sprintf("[INFO] Projecting %s into training PCA space...\n", split_name))
  rna_scores_split <- project_rna_pca(seurat_split, train_transforms$rna_fit)
  rownames(rna_scores_split) <- split_cells
  cat(sprintf("  Projected using %d genes\n", nrow(train_transforms$rna_fit$rotation)))
  
  # Build k-NN smoothing matrix WITHIN this split only
  cat(sprintf("[INFO] Building k-NN within %s cells (k=%d)...\n", split_name, K_NEIGHBORS))
  A_smooth_split <- build_knn_smooth(rna_scores_split, split_cells, 
                                      k_smooth = K_NEIGHBORS, seed = SEED_METACELL)
  
  neighbors_per_cell <- Matrix::rowSums(A_smooth_split > 0)
  cat(sprintf("  k-NN stats: mean=%.1f, min=%d, max=%d neighbors\n",
              mean(neighbors_per_cell), min(neighbors_per_cell), max(neighbors_per_cell)))
  
  # Smooth expression
  cat(sprintf("[INFO] Smoothing %s expression...\n", split_name))
  rna_counts_split <- GetAssayData(seurat_split, assay = "RNA", layer = "counts")
  smoothed_split <- smooth_rna_expression(rna_counts_split, A_smooth_split)
  cat(sprintf("[OK] %s: %d cells smoothed\n", split_name, ncol(smoothed_split$rna_log1p)))
  
  # Save smoothed data
  split_output <- list(
    rna_log1p = smoothed_split$rna_log1p,
    rna_lib = smoothed_split$rna_lib,
    rna_counts = smoothed_split$rna_counts,
    rna_pca = rna_scores_split,
    cells = split_cells,
    split_name = split_name,
    processing_seed = SEED_METACELL,
    k_smooth = K_NEIGHBORS,
    pca_dims = PCA_DIMS
  )
  saveRDS(split_output, file.path(OUTPUT_METACELLS_DIR, paste0("smoothed_", split_name, ".rds")))
  cat(sprintf("[OK] Saved smoothed_%s.rds\n", split_name))
  
  # Store for plots
  all_pca_scores[[split_name]] <- rna_scores_split
  all_splits_info[[split_name]] <- list(n_cells = length(split_cells))
}

# ============================================================================
# GENERATE PLOTS
# ============================================================================
cat("\n=== Generating plots ===\n\n")

# ----------------------------------------------------------------------------
# 1. PCA VARIANCE EXPLAINED (Scree Plot)
# ----------------------------------------------------------------------------
var_explained <- (rna_fit$sdev^2) / sum(rna_fit$sdev^2) * 100
cum_var <- cumsum(var_explained)

scree_df <- data.frame(
  PC = 1:length(var_explained),
  Variance = var_explained,
  Cumulative = cum_var
)

p_scree <- ggplot(scree_df, aes(x = PC)) +
  geom_bar(aes(y = Variance), stat = "identity", fill = colorblind_palette[2], 
           color = "black", alpha = 0.7, width = 0.7) +
  geom_line(aes(y = Cumulative/2), color = colorblind_palette[6], linewidth = 1.2) +
  geom_point(aes(y = Cumulative/2), color = colorblind_palette[6], size = 2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_y_continuous(
    name = "Variance Explained (%)",
    sec.axis = sec_axis(~.*2, name = "Cumulative Variance (%)")
  ) +
  labs(
    x = "Principal Component",
    title = "PCA Variance Explained (RNA)",
    subtitle = sprintf("%s | %d HVGs | %d PCs", SAMPLE_NAME, N_VARIABLE_FEATURES, PCA_DIMS),
    caption = sprintf("Total variance captured: %.1f%%", cum_var[PCA_DIMS])
  ) +
  theme_pub() +
  theme(
    axis.title.y.right = element_text(color = colorblind_palette[6]),
    axis.text.y.right = element_text(color = colorblind_palette[6])
  )

save_plot(p_scree, file.path(plots_dir, paste0(SAMPLE_NAME, "_pca_scree")),
                       width = 10, height = 6)

# ----------------------------------------------------------------------------
# 2. PCA EMBEDDING BY SPLIT (PC1 vs PC2)
# ----------------------------------------------------------------------------
pca_plot_df <- do.call(rbind, lapply(names(all_pca_scores), function(split) {
  scores <- all_pca_scores[[split]]
  data.frame(
    PC1 = scores[, 1],
    PC2 = scores[, 2],
    PC3 = if(ncol(scores) >= 3) scores[, 3] else NA,
    split = split,
    cell = rownames(scores)
  )
}))
pca_plot_df$split <- factor(pca_plot_df$split, levels = c("train", "validation", "test"))

# PC1 vs PC2
p_pca12 <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = split)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = split_colors_pub, name = "Data Split",
                     labels = c("Train", "Validation", "Test")) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA Embedding by Data Split",
    subtitle = sprintf("%s | RNA-only | k=%d neighbors", SAMPLE_NAME, K_NEIGHBORS)
  ) +
  theme_pub() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_plot(p_pca12, file.path(plots_dir, paste0(SAMPLE_NAME, "_pca_pc1_pc2_by_split")),
                       width = 10, height = 8)

# PC1 vs PC3
if (!all(is.na(pca_plot_df$PC3))) {
  p_pca13 <- ggplot(pca_plot_df, aes(x = PC1, y = PC3, color = split)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_color_manual(values = split_colors_pub, name = "Data Split",
                       labels = c("Train", "Validation", "Test")) +
    labs(
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC3 (%.1f%%)", var_explained[3]),
      title = "PCA Embedding (PC1 vs PC3)",
      subtitle = sprintf("%s | RNA-only", SAMPLE_NAME)
    ) +
    theme_pub() +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  save_plot(p_pca13, file.path(plots_dir, paste0(SAMPLE_NAME, "_pca_pc1_pc3_by_split")),
                         width = 10, height = 8)
}

# ----------------------------------------------------------------------------
# 3. SPLIT-FACETED PCA PLOT
# ----------------------------------------------------------------------------
p_pca_facet <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = split)) +
  geom_point(alpha = 0.5, size = 0.6) +
  scale_color_manual(values = split_colors_pub, name = "Data Split",
                     labels = c("Train", "Validation", "Test")) +
  facet_wrap(~ split, labeller = labeller(split = c(
    "train" = sprintf("Train (n=%d)", all_splits_info$train$n_cells),
    "validation" = sprintf("Validation (n=%d)", all_splits_info$validation$n_cells),
    "test" = sprintf("Test (n=%d)", all_splits_info$test$n_cells)
  ))) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA Embedding by Data Split (Faceted)",
    subtitle = sprintf("%s | RNA PCA | k=%d neighbors", SAMPLE_NAME, K_NEIGHBORS)
  ) +
  theme_pub() +
  theme(legend.position = "none")

save_plot(p_pca_facet, file.path(plots_dir, paste0(SAMPLE_NAME, "_pca_faceted_by_split")),
                       width = 14, height = 5)

# ----------------------------------------------------------------------------
# 4. DENSITY PLOTS FOR PC1 AND PC2 BY SPLIT
# ----------------------------------------------------------------------------
p_density_pc1 <- ggplot(pca_plot_df, aes(x = PC1, fill = split, color = split)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  scale_fill_manual(values = split_colors_pub, name = "Split") +
  scale_color_manual(values = split_colors_pub, name = "Split") +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = "Density",
    title = "PC1 Distribution by Split"
  ) +
  theme_pub()

p_density_pc2 <- ggplot(pca_plot_df, aes(x = PC2, fill = split, color = split)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  scale_fill_manual(values = split_colors_pub, name = "Split") +
  scale_color_manual(values = split_colors_pub, name = "Split") +
  labs(
    x = sprintf("PC2 (%.1f%%)", var_explained[2]),
    y = "Density",
    title = "PC2 Distribution by Split"
  ) +
  theme_pub()

p_density_combined <- p_density_pc1 + p_density_pc2 +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "PCA Component Distributions by Data Split",
    subtitle = sprintf("%s | Validation of split randomization", SAMPLE_NAME),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40")
    )
  )

save_plot(p_density_combined, 
                       file.path(plots_dir, paste0(SAMPLE_NAME, "_pca_density_by_split")),
                       width = 14, height = 6)

# ----------------------------------------------------------------------------
# 5. SMOOTHING EFFECT VISUALIZATION
# ----------------------------------------------------------------------------
# Compare variance before and after smoothing for a sample of genes

# Get raw and smoothed data for training
raw_train <- as.matrix(rna_counts_train[1:min(100, nrow(rna_counts_train)), ])
smoothed_train_mat <- as.matrix(smoothed_train$rna_counts[1:min(100, nrow(smoothed_train$rna_counts)), ])

# Calculate variance per gene
var_raw <- apply(raw_train, 1, var)
var_smooth <- apply(smoothed_train_mat, 1, var)

var_df <- data.frame(
  gene = names(var_raw),
  raw = var_raw,
  smoothed = var_smooth
) %>%
  mutate(var_ratio = smoothed / (raw + 1e-10))

p_var_reduction <- ggplot(var_df, aes(x = log10(raw + 1), y = log10(smoothed + 1))) +
  geom_point(alpha = 0.5, color = colorblind_palette[2], size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "Log10 Variance (Raw Counts + 1)",
    y = "Log10 Variance (Smoothed Counts + 1)",
    title = "k-NN Smoothing Effect on Gene Variance",
    subtitle = sprintf("%s | k=%d neighbors | First 100 genes", SAMPLE_NAME, K_NEIGHBORS),
    caption = sprintf("Points below diagonal = variance reduction\nMedian reduction: %.1f%%",
                      (1 - median(var_df$var_ratio, na.rm = TRUE)) * 100)
  ) +
  theme_pub() +
  coord_fixed()

save_plot(p_var_reduction, 
                       file.path(plots_dir, paste0(SAMPLE_NAME, "_smoothing_variance_reduction")),
                       width = 8, height = 8)

# ----------------------------------------------------------------------------
# 6. SUMMARY STATISTICS TABLE PLOT
# ----------------------------------------------------------------------------
summary_df <- data.frame(
  Split = c("Train", "Validation", "Test"),
  Cells = c(all_splits_info$train$n_cells, 
            all_splits_info$validation$n_cells, 
            all_splits_info$test$n_cells),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Percentage = sprintf("%.1f%%", Cells / sum(Cells) * 100),
    k_Neighbors = K_NEIGHBORS
  )

p_summary <- ggplot(summary_df, aes(x = Split, y = Cells, fill = Split)) +
  geom_col(color = "black", alpha = 0.8, width = 0.6) +
  geom_text(aes(label = format(Cells, big.mark = ",")), 
            vjust = -0.5, size = 4, fontface = "bold") +
  geom_text(aes(label = Percentage), vjust = 1.5, size = 3.5, color = "white") +
  scale_fill_manual(values = split_colors_pub) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "",
    y = "Number of Cells",
    title = "Data Split Distribution",
    subtitle = sprintf("%s | RNA-only Metacell Creation", SAMPLE_NAME),
    caption = sprintf("Method: PCA + k-NN smoothing (k=%d, %d dims)", K_NEIGHBORS, PCA_DIMS)
  ) +
  theme_pub() +
  theme(legend.position = "none")

save_plot(p_summary, 
                       file.path(plots_dir, paste0(SAMPLE_NAME, "_split_summary")),
                       width = 8, height = 6)

# ============================================================================
# COMPLETION SUMMARY
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 050 COMPLETE (scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("k-NN smoothing completed for all data splits (RNA-ONLY)\n\n")
cat("Parameters:\n")
cat(sprintf("  - k = %d neighbors\n", K_NEIGHBORS))
cat(sprintf("  - PCA dims = %d\n", PCA_DIMS))
cat(sprintf("  - Variable genes = %d\n", N_VARIABLE_FEATURES))
cat(sprintf("  - Random seed = %d\n", SEED_METACELL))
cat(sprintf("  - Method: %s\n", DIMRED_METHOD_SUFFIX))

cat("\nKey features:\n")
cat("  ✓ Training: Built PCA transformation from scratch\n")
cat("  ✓ Val/Test: Projected into training coordinate system\n")
cat("  ✓ Val/Test: Smoothed only among val/test cells (no data leakage)\n")
cat("  ✓ All cells preserved (1-to-1 mapping)\n")
cat("  ✓ RNA-ONLY: No ATAC data used\n")

cat("\nFiles saved:\n")
cat(sprintf("  - %s/train_transforms.rds\n", basename(OUTPUT_METACELLS_DIR)))
cat(sprintf("  - %s/smoothed_train.rds\n", basename(OUTPUT_METACELLS_DIR)))
cat(sprintf("  - %s/smoothed_validation.rds\n", basename(OUTPUT_METACELLS_DIR)))
cat(sprintf("  - %s/smoothed_test.rds\n", basename(OUTPUT_METACELLS_DIR)))

cat("\nPublication-ready plots saved:\n")
cat(sprintf("  - %s_pca_scree.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_pca_pc1_pc2_by_split.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_pca_pc1_pc3_by_split.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_pca_faceted_by_split.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_pca_density_by_split.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_smoothing_variance_reduction.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_split_summary.pdf/png\n", SAMPLE_NAME))

cat("\nNext: Run 04_feature_extraction.R\n")
