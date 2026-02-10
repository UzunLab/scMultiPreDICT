# ============================================================================
# Step_040: Gene-Specific Feature Extraction (scRNA-only, Automated)
# ============================================================================
# This script extracts HVG expression features for predicting target gene
# expression. Uses RNA data only (no ATAC peaks).
#
# Input: Smoothed metacell data from Step_030 (k-NN smoothing)
# Output: Gene-specific feature matrices for train/val/test splits
#
# Features:
#   - HVG expression as predictor features
#   - Target gene excluded from its own features (prevents leakage)
#   - HVGs computed from TRAINING data only
#   - Publication-ready diagnostic plots
#
# Usage:
#   1. Edit config.R with your parameters
#   2. Run: Rscript 04_feature_extraction.R
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
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(Cairo)
  library(viridis)
})

# If publication helpers are not available, define safe fallbacks here.
if (!exists("theme_pub") || !is.function(theme_pub)) {
  theme_pub <- function(base_size = 12) {
    theme_minimal(base_size = base_size) +
      theme(
        plot.title = element_text(face = "bold", size = base_size + 2),
        plot.subtitle = element_text(size = base_size),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "right"
      )
  }
}

if (!exists("save_plot") || !is.function(save_plot)) {
  save_plot <- function(plot, filename, width = 7, height = 5, dpi = 600) {
    pdf_file <- paste0(filename, ".pdf")
    png_file <- paste0(filename, ".png")
    # Try to save PDF via Cairo where available
    try({
      cairo_pdf(pdf_file, width = width, height = height)
      print(plot)
      dev.off()
    }, silent = TRUE)
    # Fallback to ggsave PNG
    try({
      ggsave(png_file, plot = plot, width = width, height = height, dpi = dpi)
    }, silent = TRUE)
    cat(sprintf("Saved: %s, %s\n", pdf_file, png_file))
  }
}

# ============================================================================
# PLOTTING SETUP
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

# Gene set variables and loading logic are handled later after data is loaded.

# ============================================================================
# PRINT CONFIGURATION
# ============================================================================
print_config()
print_output_directories()

# Create output directories
create_output_directories()

# Create plots directory
plots_dir <- file.path(OUTPUT_FEATURES_DIR, "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 060: Feature Extraction (scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat(sprintf("Dimensionality Reduction Method: %s\n", DIM_REDUCTION_METHOD))
cat(sprintf("Output directory: %s\n\n", OUTPUT_FEATURES_DIR))

set.seed(SEED_FEATURES)

# ============================================================================
# LOAD INPUT DATA
# ============================================================================
cat("=== Loading preprocessed data ===\n\n")

# Load the original Seurat object (for HVG selection)
seurat_obj <- readRDS(path.expand(INPUT_SEURAT_SPLITS))
cat(sprintf("Loaded Seurat object: %d cells, %d genes\n",
            ncol(seurat_obj), nrow(seurat_obj[["RNA"]])))

# Load smoothed RNA data for all splits
cat("\nLoading smoothed RNA data from metacell creation...\n")
smoothed_train <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_train.rds"))
smoothed_val <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_validation.rds"))
smoothed_test <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_test.rds"))

cat(sprintf("Train: %d cells\n", length(smoothed_train$cells)))
cat(sprintf("Validation: %d cells\n", length(smoothed_val$cells)))
cat(sprintf("Test: %d cells\n", length(smoothed_test$cells)))

# ============================================================================

# Extract training cells only (always needed for fallback HVG computation)
train_cells_for_hvg <- colnames(seurat_obj)[seurat_obj$data_split == "train"]
seurat_train_subset <- seurat_obj[, train_cells_for_hvg]

cat(sprintf("Using %d training cells (excluding %d val + %d test cells)\n",
            length(train_cells_for_hvg),
            sum(seurat_obj$data_split == "validation"),
            sum(seurat_obj$data_split == "test")))

all_hvgs <- NULL
top_hvg_features <- NULL
hvg_genes <- NULL
random_genes <- NULL

# Interpret user-provided lists as TARGET gene lists. Predictor features (top HVGs)
# are computed from TRAINING data (no leakage) using FindVariableFeatures.
hvg_present <- exists("HVG_GENE_FILE") && !is.null(HVG_GENE_FILE) && HVG_GENE_FILE != "" && file.exists(path.expand(HVG_GENE_FILE))
rand_present <- exists("RANDOM_GENE_FILE") && !is.null(RANDOM_GENE_FILE) && RANDOM_GENE_FILE != "" && file.exists(path.expand(RANDOM_GENE_FILE))

if (!hvg_present && !rand_present) {
  stop("No target gene list provided. Please set HVG_GENE_FILE and/or RANDOM_GENE_FILE in config.R (these are treated as TARGET gene lists).")
}

# Compute top HVG features from training Seurat subset (training-only, avoids leakage)
DefaultAssay(seurat_train_subset) <- "RNA"
cat("Computing top HVG features from training data (no leakage)...\n")
seurat_train_subset <- FindVariableFeatures(
  seurat_train_subset,
  selection.method = "vst",
  nfeatures = N_HVG_FEATURES,
  verbose = FALSE
)
all_hvgs <- VariableFeatures(seurat_train_subset)
top_hvg_features <- head(all_hvgs, min(N_HVG_FEATURES, length(all_hvgs)))
cat(sprintf("✓ Computed %d top HVG predictor features from training data\n", length(top_hvg_features)))

# Load target gene lists (user-supplied)
if (hvg_present) {
  tgt <- readLines(path.expand(HVG_GENE_FILE))
  tgt <- tgt[nchar(tgt) > 0]
  hvg_genes <- unique(tgt)
  cat(sprintf("Loaded %d target genes from HVG_GENE_FILE (treated as targets): %s\n", length(hvg_genes), basename(HVG_GENE_FILE)))
}

if (rand_present) {
  random_genes <- readLines(path.expand(RANDOM_GENE_FILE))
  random_genes <- unique(random_genes[nchar(random_genes) > 0])
  cat(sprintf("Loaded %d target genes from RANDOM_GENE_FILE: %s\n", length(random_genes), basename(RANDOM_GENE_FILE)))
}

# Create list of gene sets to process (these are TARGET sets)
gene_sets <- list()
if (!is.null(random_genes) && length(random_genes) > 0) gene_sets$Random_genes <- random_genes
if (!is.null(hvg_genes) && length(hvg_genes) > 0) gene_sets$HVG <- hvg_genes

cat(sprintf("\nTotal gene sets to process: %d\n", length(gene_sets)))
for (set_name in names(gene_sets)) {
  cat(sprintf("  - %s: %d genes\n", set_name, length(gene_sets[[set_name]])))
}

# Assign colors for each gene set for plotting
geneset_colors <- setNames(colorblind_palette[seq_len(length(gene_sets))], names(gene_sets))

# ============================================================================
# FEATURE EXTRACTION FUNCTIONS
# ============================================================================

# Extract HVG expression features (excluding target gene)
extract_hvg_features <- function(gene_name, hvg_list, smoothed_data) {
  hvg_features <- setdiff(hvg_list, gene_name)
  all_gene_names <- rownames(smoothed_data$rna_log1p)
  hvg_features <- intersect(hvg_features, all_gene_names)
  
  if (length(hvg_features) == 0) {
    warning(sprintf("No HVG features available for gene %s", gene_name))
    return(NULL)
  }
  
  hvg_matrix <- smoothed_data$rna_log1p[hvg_features, , drop = FALSE]
  hvg_matrix <- t(hvg_matrix)  # cells × genes
  return(hvg_matrix)
}

# Extract target gene expression
extract_gene_expression <- function(gene_name, smoothed_data) {
  all_gene_names <- rownames(smoothed_data$rna_log1p)
  
  if (!gene_name %in% all_gene_names) {
    warning(sprintf("Gene %s not found in RNA data", gene_name))
    return(NULL)
  }
  
  expression <- smoothed_data$rna_log1p[gene_name, ]
  as.numeric(expression)
}

# ============================================================================
# PROCESS EACH GENE SET
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("PROCESSING GENE SETS (RNA-ONLY)\n")
cat("=", rep("=", 70), "\n", sep = "")

# Storage for summary across gene sets
all_gene_sets_summary <- list()

for (set_name in names(gene_sets)) {
  cat(sprintf("\n\n>>> Processing gene set: %s (%d genes) <<<\n", 
              set_name, length(gene_sets[[set_name]])))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  target_genes <- gene_sets[[set_name]]
  
  # Create output directory for this gene set
  set_output_dir <- file.path(OUTPUT_FEATURES_DIR, set_name)
  if (!dir.exists(set_output_dir)) {
    dir.create(set_output_dir, recursive = TRUE)
    cat(sprintf("Created output directory: %s\n", set_output_dir))
  }
  
  # Check which target genes exist in the RNA data
  all_gene_names <- rownames(smoothed_train$rna_log1p)
  valid_target_genes <- intersect(target_genes, all_gene_names)
  missing_genes <- setdiff(target_genes, all_gene_names)
  
  if (length(missing_genes) > 0) {
    cat(sprintf("\nWarning: %d target genes not found in RNA data:\n", length(missing_genes)))
    cat(paste(head(missing_genes, 10), collapse = ", "), "\n")
    if (length(missing_genes) > 10) cat("... and", length(missing_genes) - 10, "more\n")
  }
  
  if (length(valid_target_genes) == 0) {
    cat(sprintf("ERROR: No genes in %s set found in RNA data. Skipping.\n", set_name))
    next
  }
  
  cat(sprintf("\nTarget genes found in RNA data: %d\n", length(valid_target_genes)))
  
  # Extract features for all genes
  cat("\n=== Extracting HVG features for each gene ===\n")
  
  gene_features <- list()
  n_genes <- length(valid_target_genes)
  progress_interval <- max(1, floor(n_genes / 20))
  
  for (i in seq_along(valid_target_genes)) {
    gene_name <- valid_target_genes[i]
    
    if (i %% progress_interval == 0 || i == n_genes) {
      cat(sprintf("  [%d/%d] %.1f%% complete\r", i, n_genes, 100 * i / n_genes))
    }
    
    tryCatch({
      # Extract features
      train_hvgs <- extract_hvg_features(gene_name, top_hvg_features, smoothed_train)
      val_hvgs <- extract_hvg_features(gene_name, top_hvg_features, smoothed_val)
      test_hvgs <- extract_hvg_features(gene_name, top_hvg_features, smoothed_test)
      
      # Extract targets
      train_y <- extract_gene_expression(gene_name, smoothed_train)
      val_y <- extract_gene_expression(gene_name, smoothed_val)
      test_y <- extract_gene_expression(gene_name, smoothed_test)
      
      gene_features[[gene_name]] <- list(
        gene_name = gene_name,
        n_hvg_features = ncol(train_hvgs),
        n_total_features = ncol(train_hvgs),
        hvg_feature_names = colnames(train_hvgs),
        
        train = list(
          X = train_hvgs,
          y = train_y,
          cells = smoothed_train$cells
        ),
        validation = list(
          X = val_hvgs,
          y = val_y,
          cells = smoothed_val$cells
        ),
        test = list(
          X = test_hvgs,
          y = test_y,
          cells = smoothed_test$cells
        )
      )
    }, error = function(e) {
      cat(sprintf("\nWarning: Failed for gene %s: %s\n", gene_name, e$message))
    })
  }
  
  cat("\n")
  
  # Remove any failed genes
  gene_features <- gene_features[!sapply(gene_features, is.null)]
  
  # Summary
  cat(sprintf("\n=== Feature Extraction Summary for %s ===\n", set_name))
  cat(sprintf("Successfully extracted features for %d genes\n", length(gene_features)))
  
  if (length(gene_features) > 0) {
    n_hvg_features_vec <- sapply(gene_features, function(x) x$n_hvg_features)
    
    cat(sprintf("\nHVG features per gene:\n"))
    cat(sprintf("  Mean: %.1f\n", mean(n_hvg_features_vec)))
    cat(sprintf("  Range: [%d, %d]\n", min(n_hvg_features_vec), max(n_hvg_features_vec)))
    
    # Example
    first_gene <- gene_features[[1]]
    cat(sprintf("\nExample: Gene '%s'\n", first_gene$gene_name))
    cat(sprintf("  Train: %d cells × %d features\n", 
                nrow(first_gene$train$X), ncol(first_gene$train$X)))
  }
  
  # Save features
  cat("\n=== Saving extracted features ===\n")
  
  output_file <- file.path(set_output_dir, "gene_specific_features.rds")
  saveRDS(gene_features, output_file)
  cat(sprintf("Saved: %s\n", output_file))
  
  # Save metadata
  if (length(gene_features) > 0) {
    feature_metadata <- data.frame(
      gene_name = names(gene_features),
      n_hvg_features = sapply(gene_features, function(x) x$n_hvg_features),
      n_total_features = sapply(gene_features, function(x) x$n_total_features),
      n_train_cells = sapply(gene_features, function(x) nrow(x$train$X)),
      n_val_cells = sapply(gene_features, function(x) nrow(x$validation$X)),
      n_test_cells = sapply(gene_features, function(x) nrow(x$test$X)),
      mean_train_expr = sapply(gene_features, function(x) mean(x$train$y)),
      var_train_expr = sapply(gene_features, function(x) var(x$train$y)),
      stringsAsFactors = FALSE
    )
    
    metadata_file <- file.path(set_output_dir, "gene_features_metadata.csv")
    write.csv(feature_metadata, metadata_file, row.names = FALSE)
    cat(sprintf("Saved metadata: %s\n", metadata_file))
    
    all_gene_sets_summary[[set_name]] <- feature_metadata
  }
  
  # Save parameters
  params <- list(
    gene_set = set_name,
    feature_type = "RNA_only",
    n_hvg_features = N_HVG_FEATURES,
    n_genes_extracted = length(gene_features),
    seed = SEED_FEATURES,
    extraction_date = Sys.time(),
    dim_reduction_method = DIM_REDUCTION_METHOD
  )
  
  params_file <- file.path(set_output_dir, "feature_extraction_params.rds")
  saveRDS(params, params_file)
  cat(sprintf("Saved parameters: %s\n", params_file))
  
  cat(sprintf("\n>>> Completed %s gene set <<<\n", set_name))
}

# ============================================================================
# GENERATE PLOTS
# ============================================================================
cat("\n=== Generating plots ===\n\n")

# Combine all gene set summaries
if (length(all_gene_sets_summary) > 0) {
  
  combined_summary <- do.call(rbind, lapply(names(all_gene_sets_summary), function(set) {
    df <- all_gene_sets_summary[[set]]
    df$gene_set <- set
    df
  }))
  
  # ----------------------------------------------------------------------------
  # 1. TARGET GENE EXPRESSION DISTRIBUTION
  # ----------------------------------------------------------------------------
  p_expr_dist <- ggplot(combined_summary, aes(x = mean_train_expr, fill = gene_set)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "black", position = "identity") +
    scale_fill_manual(values = geneset_colors, name = "Gene Set") +
    labs(
      x = "Mean Expression (log1p CPM)",
      y = "Number of Genes",
      title = "Target Gene Expression Distribution",
      subtitle = sprintf("%s | Training data", SAMPLE_NAME)
    ) +
    theme_pub() +
    facet_wrap(~ gene_set, ncol = 2)
  
  save_plot(p_expr_dist, 
            file.path(plots_dir, paste0(SAMPLE_NAME, "_target_expression_distribution")),
            width = 12, height = 5)
  
  # ----------------------------------------------------------------------------
  # 2. EXPRESSION VARIANCE DISTRIBUTION
  # ----------------------------------------------------------------------------
  p_var_dist <- ggplot(combined_summary, aes(x = log10(var_train_expr + 0.01), fill = gene_set)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "black", position = "identity") +
    scale_fill_manual(values = geneset_colors, name = "Gene Set") +
    labs(
      x = "Log10 Expression Variance",
      y = "Number of Genes",
      title = "Target Gene Variance Distribution",
      subtitle = sprintf("%s | Training data", SAMPLE_NAME)
    ) +
    theme_pub() +
    facet_wrap(~ gene_set, ncol = 2)
  
  save_plot(p_var_dist, 
            file.path(plots_dir, paste0(SAMPLE_NAME, "_target_variance_distribution")),
            width = 12, height = 5)
  
  # ----------------------------------------------------------------------------
  # 3. MEAN VS VARIANCE SCATTER
  # ----------------------------------------------------------------------------
  p_mean_var <- ggplot(combined_summary, aes(x = mean_train_expr, y = log10(var_train_expr + 0.01), 
                                             color = gene_set)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = geneset_colors, name = "Gene Set") +
    labs(
      x = "Mean Expression (log1p CPM)",
      y = "Log10 Expression Variance",
      title = "Target Gene Mean vs Variance",
      subtitle = sprintf("%s | Training data", SAMPLE_NAME)
    ) +
    theme_pub()
  
  save_plot(p_mean_var, 
            file.path(plots_dir, paste0(SAMPLE_NAME, "_target_mean_vs_variance")),
            width = 10, height = 8)
  
  # ----------------------------------------------------------------------------
  # 4. FEATURE COUNT SUMMARY
  # ----------------------------------------------------------------------------
  feature_summary <- combined_summary %>%
    group_by(gene_set) %>%
    summarise(
      n_genes = n(),
      mean_features = mean(n_hvg_features),
      min_features = min(n_hvg_features),
      max_features = max(n_hvg_features),
      .groups = "drop"
    )
  
  p_feature_count <- ggplot(feature_summary, aes(x = gene_set, y = n_genes, fill = gene_set)) +
    geom_col(color = "black", alpha = 0.8, width = 0.6) +
    geom_text(aes(label = n_genes), vjust = -0.5, size = 5, fontface = "bold") +
    scale_fill_manual(values = geneset_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      x = "Gene Set",
      y = "Number of Target Genes",
      title = "Target Genes per Gene Set",
      subtitle = sprintf("%s | HVG features: %d", SAMPLE_NAME, N_HVG_FEATURES)
    ) +
    theme_pub() +
    theme(legend.position = "none")
  
  save_plot(p_feature_count, 
            file.path(plots_dir, paste0(SAMPLE_NAME, "_target_genes_per_set")),
            width = 8, height = 6)
  
  # ----------------------------------------------------------------------------
  # 5. SAMPLE FEATURE MATRIX HEATMAP (first 20 genes, first 50 cells)
  # ----------------------------------------------------------------------------
  # Get first gene set's features
  first_set <- names(all_gene_sets_summary)[1]
  features_file <- file.path(OUTPUT_FEATURES_DIR, first_set, "gene_specific_features.rds")
  
  if (file.exists(features_file)) {
    sample_features <- readRDS(features_file)
    
    if (length(sample_features) > 0) {
      # Get first gene's feature matrix
      first_gene <- sample_features[[1]]
      sample_X <- first_gene$train$X[1:min(50, nrow(first_gene$train$X)), 
                                     1:min(20, ncol(first_gene$train$X))]
      
      # Convert to long format for heatmap
      heatmap_df <- as.data.frame(sample_X) %>%
        mutate(cell = row_number()) %>%
        tidyr::pivot_longer(-cell, names_to = "feature", values_to = "expression")
      
      p_heatmap <- ggplot(heatmap_df, aes(x = feature, y = cell, fill = expression)) +
        geom_tile() +
        scale_fill_viridis_c(option = "viridis", name = "Expression") +
        labs(
          x = "HVG Feature",
          y = "Cell",
          title = sprintf("Sample Feature Matrix: %s", first_gene$gene_name),
          subtitle = sprintf("First 50 cells × 20 HVG features")
        ) +
        theme_pub() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      
      save_plot(p_heatmap, 
                file.path(plots_dir, paste0(SAMPLE_NAME, "_sample_feature_matrix")),
                width = 12, height = 8)
    }
  }
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 040 COMPLETE (scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Gene-specific HVG expression features extracted successfully\n\n")

cat("Summary by gene set:\n")
for (set_name in names(gene_sets)) {
  set_output_dir <- file.path(OUTPUT_FEATURES_DIR, set_name)
  feature_file <- file.path(set_output_dir, "gene_specific_features.rds")
  
  if (file.exists(feature_file)) {
    features <- readRDS(feature_file)
    cat(sprintf("  %s: %d genes processed\n", set_name, length(features)))
  }
}

cat("\nKey parameters:\n")
cat(sprintf("  - Feature type: RNA expression only (HVGs)\n"))
cat(sprintf("  - Top HVG features: %d\n", N_HVG_FEATURES))
cat(sprintf("  - Dimensionality reduction: %s\n", DIMRED_METHOD_SUFFIX))
cat(sprintf("  - Random seed: %d\n", SEED_FEATURES))
cat(sprintf("  - Total train cells: %d\n", length(smoothed_train$cells)))
cat(sprintf("  - Total validation cells: %d\n", length(smoothed_val$cells)))
cat(sprintf("  - Total test cells: %d\n", length(smoothed_test$cells)))

cat("\nOutput structure:\n")
for (set_name in names(gene_sets)) {
  cat(sprintf("  %s/\n", set_name))
  cat("    - gene_specific_features.rds\n")
  cat("    - gene_features_metadata.csv\n")
  cat("    - feature_extraction_params.rds\n")
}

cat("\nPublication-ready plots saved:\n")
cat(sprintf("  - %s_target_expression_distribution.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_target_variance_distribution.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_target_mean_vs_variance.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_target_genes_per_set.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_sample_feature_matrix.pdf/png\n", SAMPLE_NAME))

cat("\nIMPORTANT NOTES:\n")
cat("  ✓ Target gene is EXCLUDED from HVG features to prevent data leakage\n")
cat("  ✓ HVGs computed from TRAINING data only\n")
cat("  ✓ All feature matrices are in (cells × features) format\n")

cat("\nNext: Run 05_linear_tree_models.R\n")
