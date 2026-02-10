# ============================================================================
# Step_02b: Select Target Genes for Prediction (OPTIONAL - Standalone Script)
# ============================================================================
# Purpose: Auto-select random non-HVG genes to use as prediction targets.
#
# *** This script is NOT part of the default pipeline. ***
# It is only needed when analyzing NEW datasets that do not have
# pre-computed target gene lists in data/target_genes/.
# For published datasets, use the pre-computed files instead.
#
# Before running this script, uncomment and set these parameters in config.R:
#   N_RANDOM_TARGET_GENES <- 100
#   TARGET_MIN_DETECTION  <- 5
#   TARGET_MAX_DETECTION  <- 95
#   TARGET_MIN_MEAN_EXPR  <- 0.1
#   SEED_TARGET_GENES     <- 2025
#   OUTPUT_TARGET_GENES_DIR <- file.path(BASE_OUTPUT_DIR, "target_genes", SAMPLE_NAME)
#
# Strategy:
# - Exclude highly variable genes (HVGs) used for PCA/dimensionality reduction
# - Select genes with reasonable expression levels
# - Avoid very rare or very sparse genes
# - Save gene list for downstream prediction tasks
#
# Input: Seurat object with split labels from Step 02a
# Output: Target gene list (text file) and summary (CSV)
#
# Usage:
#   cd combined
#   Rscript R/02b_select_target_genes.R
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
  library(Matrix)
  library(ggplot2)
})

# Create output directories
create_output_directories()

# Set seed for reproducibility
set.seed(SEED_TARGET_GENES)

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 035: Select Target Genes for Prediction\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Parameters:\n")
cat(sprintf("  Number of random target genes: %d\n", N_RANDOM_TARGET_GENES))
cat(sprintf("  Detection rate range: %.0f%% - %.0f%%\n", TARGET_MIN_DETECTION, TARGET_MAX_DETECTION))
cat(sprintf("  Minimum mean expression: %.2f\n", TARGET_MIN_MEAN_EXPR))
cat(sprintf("  Random seed: %d\n\n", SEED_TARGET_GENES))

# ============================================================================
# LOAD DATA
# ============================================================================
cat("=== Loading Seurat object ===\n\n")

input_file <- file.path(
  OUTPUT_SPLITS_DIR,
  paste0(SAMPLE_NAME, "_seurat_obj_with_splits.rds")
)

if (!file.exists(input_file)) {
  stop(sprintf("ERROR: Input file not found: %s\nPlease run Step_020 first.", input_file))
}

cat(sprintf("Loading: %s\n", input_file))
seurat_obj <- readRDS(input_file)

cat(sprintf("Total cells: %d\n", ncol(seurat_obj)))
cat(sprintf("Total genes: %d\n", nrow(seurat_obj[["RNA"]])))

# ============================================================================
# GET TRAINING CELLS ONLY
# ============================================================================
# We compute statistics on training data only to avoid data leakage
train_cells <- colnames(seurat_obj)[seurat_obj$data_split == "train"]
seurat_train <- seurat_obj[, train_cells]

cat(sprintf("Training cells: %d\n\n", length(train_cells)))

# ============================================================================
# FIND HIGHLY VARIABLE GENES (to exclude)
# ============================================================================
cat("=== Finding highly variable genes ===\n\n")

DefaultAssay(seurat_train) <- "RNA"
seurat_train <- NormalizeData(seurat_train, verbose = FALSE)
seurat_train <- FindVariableFeatures(seurat_train, nfeatures = N_VARIABLE_FEATURES, verbose = FALSE)

hvgs <- VariableFeatures(seurat_train)
cat(sprintf("HVGs (excluded from random selection): %d\n", length(hvgs)))

# Get all genes
all_genes <- rownames(seurat_train[["RNA"]])

# Exclude HVGs
non_hvgs <- setdiff(all_genes, hvgs)
cat(sprintf("Non-HVG genes available: %d\n\n", length(non_hvgs)))

# ============================================================================
# CALCULATE GENE STATISTICS
# ============================================================================
cat("=== Calculating gene statistics ===\n\n")

rna_counts <- GetAssayData(seurat_train, assay = "RNA", layer = "counts")
rna_counts_non_hvg <- rna_counts[non_hvgs, , drop = FALSE]

# Compute detection rate (% cells expressing the gene) and mean expression
detection_rate <- Matrix::rowMeans(rna_counts_non_hvg > 0) * 100
mean_expression <- Matrix::rowMeans(rna_counts_non_hvg)

# Create statistics dataframe for all non-HVG genes
gene_stats <- data.frame(
  gene = non_hvgs,
  detection_rate = detection_rate,
  mean_expression = mean_expression,
  stringsAsFactors = FALSE
)

cat(sprintf("Gene statistics computed for %d non-HVG genes\n", nrow(gene_stats)))
cat(sprintf("  Detection rate - range: %.1f%% - %.1f%%\n", 
            min(gene_stats$detection_rate), max(gene_stats$detection_rate)))
cat(sprintf("  Mean expression - range: %.4f - %.2f\n\n",
            min(gene_stats$mean_expression), max(gene_stats$mean_expression)))

# ============================================================================
# FILTER ELIGIBLE GENES
# ============================================================================
cat("=== Filtering eligible genes ===\n\n")

# Apply filters
eligible_mask <- gene_stats$detection_rate >= TARGET_MIN_DETECTION &
                 gene_stats$detection_rate <= TARGET_MAX_DETECTION &
                 gene_stats$mean_expression >= TARGET_MIN_MEAN_EXPR

eligible_genes <- gene_stats$gene[eligible_mask]

cat(sprintf("Eligibility criteria:\n"))
cat(sprintf("  Detection rate: %.0f%% - %.0f%%\n", TARGET_MIN_DETECTION, TARGET_MAX_DETECTION))
cat(sprintf("  Mean expression: > %.2f\n", TARGET_MIN_MEAN_EXPR))
cat(sprintf("Eligible genes: %d / %d (%.1f%%)\n\n", 
            length(eligible_genes), length(non_hvgs),
            100 * length(eligible_genes) / length(non_hvgs)))

# ============================================================================
# SELECT RANDOM TARGET GENES
# ============================================================================
cat("=== Selecting random target genes ===\n\n")

n_targets <- N_RANDOM_TARGET_GENES

if (length(eligible_genes) < n_targets) {
  warning(sprintf("Only %d eligible genes available. Using all of them.", length(eligible_genes)))
  target_genes <- eligible_genes
} else {
  target_genes <- sample(eligible_genes, n_targets)
}

cat(sprintf("Selected %d target genes for prediction\n\n", length(target_genes)))

# Create summary data frame
target_summary <- data.frame(
  gene = target_genes,
  detection_rate = detection_rate[target_genes],
  mean_expression = mean_expression[target_genes],
  is_hvg = FALSE,
  stringsAsFactors = FALSE
)

# Sort by gene name
target_summary <- target_summary[order(target_summary$gene), ]

# ============================================================================
# DISPLAY SUMMARY
# ============================================================================
cat("=== Target gene summary ===\n\n")

cat(sprintf("Detection rate:\n"))
cat(sprintf("  Mean: %.1f%%\n", mean(target_summary$detection_rate)))
cat(sprintf("  Median: %.1f%%\n", median(target_summary$detection_rate)))
cat(sprintf("  Range: %.1f%% - %.1f%%\n\n",
            min(target_summary$detection_rate), max(target_summary$detection_rate)))

cat(sprintf("Mean expression:\n"))
cat(sprintf("  Mean: %.3f\n", mean(target_summary$mean_expression)))
cat(sprintf("  Median: %.3f\n", median(target_summary$mean_expression)))
cat(sprintf("  Range: %.3f - %.3f\n\n",
            min(target_summary$mean_expression), max(target_summary$mean_expression)))

# ============================================================================
# SAVE OUTPUTS
# ============================================================================
cat("=== Saving outputs ===\n\n")

# Create output directory if needed
output_dir <- path.expand(OUTPUT_TARGET_GENES_DIR)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define output filenames
output_file_list <- file.path(output_dir, paste0("target_genes_random_", n_targets, ".txt"))
output_file_summary <- file.path(output_dir, paste0("target_genes_random_", n_targets, "_summary.csv"))
output_file_stats <- file.path(output_dir, "all_non_hvg_gene_stats.csv")

# Save target gene list (simple text file)
writeLines(target_genes, output_file_list)
cat(sprintf("✓ Saved target gene list: %s\n", output_file_list))

# Save target gene summary
write.csv(target_summary, output_file_summary, row.names = FALSE)
cat(sprintf("✓ Saved target gene summary: %s\n", output_file_summary))

# Save all non-HVG gene statistics (useful for later reference)
write.csv(gene_stats, output_file_stats, row.names = FALSE)
cat(sprintf("✓ Saved all gene statistics: %s\n", output_file_stats))

# ============================================================================
# DIAGNOSTIC PLOTS
# ============================================================================
cat("\n=== Generating diagnostic plots ===\n\n")

plots_dir <- file.path(output_dir, "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Plot 1: Detection rate distribution
p1 <- ggplot(gene_stats, aes(x = detection_rate)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = TARGET_MIN_DETECTION, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = TARGET_MAX_DETECTION, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Detection Rate Distribution (Non-HVG Genes)",
       subtitle = sprintf("Red lines: selection range (%.0f%% - %.0f%%)", 
                         TARGET_MIN_DETECTION, TARGET_MAX_DETECTION),
       x = "Detection Rate (%)", y = "Number of Genes") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(plots_dir, "01_detection_rate_distribution.png"), p1, 
       width = 8, height = 5, dpi = 600)
cat("  ✓ Saved detection rate distribution plot\n")

# Plot 2: Mean expression distribution (log scale)
p2 <- ggplot(gene_stats, aes(x = mean_expression)) +
  geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7, color = "white") +
  geom_vline(xintercept = TARGET_MIN_MEAN_EXPR, color = "red", linetype = "dashed", linewidth = 1) +
  scale_x_log10() +
  labs(title = "Mean Expression Distribution (Non-HVG Genes)",
       subtitle = sprintf("Red line: minimum threshold (%.2f)", TARGET_MIN_MEAN_EXPR),
       x = "Mean Expression (log10)", y = "Number of Genes") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(plots_dir, "02_mean_expression_distribution.png"), p2, 
       width = 8, height = 5, dpi = 600)
cat("  ✓ Saved mean expression distribution plot\n")

# Plot 3: Detection vs Expression scatter
gene_stats$selected <- gene_stats$gene %in% target_genes

p3 <- ggplot(gene_stats, aes(x = mean_expression, y = detection_rate, color = selected)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_x_log10() +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red"),
                     labels = c("Not selected", "Selected"),
                     name = "Status") +
  geom_hline(yintercept = c(TARGET_MIN_DETECTION, TARGET_MAX_DETECTION), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_vline(xintercept = TARGET_MIN_MEAN_EXPR, 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  labs(title = "Gene Selection Space",
       subtitle = sprintf("Selected %d random target genes (red)", length(target_genes)),
       x = "Mean Expression (log10)", y = "Detection Rate (%)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave(file.path(plots_dir, "03_gene_selection_scatter.png"), p3, 
       width = 8, height = 6, dpi = 600)
cat("  ✓ Saved gene selection scatter plot\n")

# Plot 4: Selected genes statistics
p4 <- ggplot(target_summary, aes(x = detection_rate)) +
  geom_histogram(bins = 20, fill = "coral", alpha = 0.7, color = "white") +
  labs(title = "Detection Rate of Selected Target Genes",
       subtitle = sprintf("n = %d genes", nrow(target_summary)),
       x = "Detection Rate (%)", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(plots_dir, "04_selected_genes_detection.png"), p4, 
       width = 7, height = 5, dpi = 600)
cat("  ✓ Saved selected genes detection plot\n")

# Combined PDF
cat("\nCreating combined PDF report...\n")

pdf(file.path(plots_dir, "02b_Target_Genes_Report.pdf"), width = 10, height = 8)

plot.new()
text(0.5, 0.6, "Step 02b: Target Gene Selection Report", cex = 2, font = 2)
text(0.5, 0.45, sprintf("Sample: %s", SAMPLE_NAME), cex = 1.2)
text(0.5, 0.35, sprintf("Selected %d random non-HVG target genes", length(target_genes)), cex = 1.2)
text(0.5, 0.25, sprintf("Date: %s", Sys.Date()), cex = 1.2)

print(p1)
print(p2)
print(p3)
print(p4)

dev.off()
cat("  ✓ Saved combined PDF report\n")

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 035 COMPLETE\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Summary:\n")
cat(sprintf("  Total genes: %d\n", length(all_genes)))
cat(sprintf("  HVGs (excluded): %d\n", length(hvgs)))
cat(sprintf("  Non-HVG genes: %d\n", length(non_hvgs)))
cat(sprintf("  Eligible genes: %d\n", length(eligible_genes)))
cat(sprintf("  Selected target genes: %d\n\n", length(target_genes)))

cat("Output files:\n")
cat(sprintf("  - %s\n", basename(output_file_list)))
cat(sprintf("  - %s\n", basename(output_file_summary)))
cat(sprintf("  - %s\n", basename(output_file_stats)))
cat(sprintf("  - plots/ (4 diagnostic plots + PDF)\n\n"))

cat("Output directory:", output_dir, "\n\n")


cat("Next step: Run Step_030 (Metacell Creation) then Step_040 (Feature Extraction)\n")
