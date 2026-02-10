# ============================================================================
# Step_01: Data Preprocessing and Quality Control (Automated)
# ============================================================================
# This script processes a single multiome dataset (RNA + ATAC) and applies
# quality control filters. It reads configuration from config.R.
#
# Input: Raw 10X Genomics multiome data (MTX format + fragments)
# Output: Filtered Seurat object with RNA and ATAC assays
#
# Usage:
#   1. Edit config.R with your dataset parameters
#   2. Run: Rscript 01_quality_control.R
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
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(Signac)
  library(readr)
  library(Matrix)
  library(rtracklayer)
  library(GenomeInfoDb)
  library(GenomicRanges)
})

# Load species-specific annotation packages
load_annotation_packages()

# Print and validate configuration
print_config()
validate_config()

# Create output directories
create_output_directories()

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 01: Quality Control\n")
cat("=", rep("=", 70), "\n\n", sep = "")

# ============================================================================
# LOAD DATA
# ============================================================================
cat("=== Loading input data ===\n\n")

cat("Loading count matrix...\n")
cat("  MTX:", INPUT_MTX, "\n")
cat("  Features:", INPUT_FEATURES, "\n")
cat("  Barcodes:", INPUT_BARCODES, "\n")

counts_all <- ReadMtx(
  mtx = INPUT_MTX,
  features = INPUT_FEATURES,
  cells = INPUT_BARCODES
)

cat(sprintf("Loaded matrix: %d features × %d cells\n", nrow(counts_all), ncol(counts_all)))

# Load features file to split RNA and ATAC
feat <- read_tsv(INPUT_FEATURES, col_names = FALSE, show_col_types = FALSE)

# Split into RNA and ATAC by feature type
# The features.tsv is expected to have at least 3 columns: feature_id, feature_name, feature_type
rna_rows <- which(feat$X3 %in% c("Gene Expression", "Gene expression", "Gene"))
atac_rows <- which(feat$X3 %in% c("Peaks", "Peak"))

cat(sprintf("RNA features: %d\n", length(rna_rows)))
cat(sprintf("ATAC features: %d\n", length(atac_rows)))

if (length(rna_rows) == 0) {
  stop("No RNA features found in features file (check column 3 values)")
}

rna_counts <- counts_all[rna_rows, , drop = FALSE]
atac_counts <- counts_all[atac_rows, , drop = FALSE]

# ============================================================================
# CREATE SEURAT OBJECT
# ============================================================================
cat("\n=== Creating Seurat object ===\n\n")

# Create Seurat Object for RNA
seurat_obj <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  project = SAMPLE_NAME
)

cat(sprintf("Created Seurat object: %d genes × %d cells\n", 
            nrow(seurat_obj), ncol(seurat_obj)))

# ============================================================================
# ADD ATAC ASSAY
# ============================================================================
cat("\n=== Adding ATAC assay ===\n\n")

# Get gene annotations
cat("Extracting gene annotations from EnsDb...\n")
ensdb <- get_ensdb()
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)

# Convert to UCSC chromosome naming style
seqlevelsStyle(annotations) <- "UCSC"
annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")
genome(annotations) <- GENOME

# Create fragment object
cat("Creating fragment object...\n")
fragments_path <- path.expand(INPUT_FRAGMENTS)
frags <- CreateFragmentObject(
  path = fragments_path, 
  cells = colnames(seurat_obj), 
  validate.fragments = TRUE
)

# Create chromatin assay and add to the seurat_obj
cat("Creating chromatin assay...\n")
seurat_obj[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  fragments = frags,
  genome = GENOME,
  sep = c(":", "-"),
  annotation = annotations
)

cat(sprintf("Added ATAC assay: %d peaks\n", nrow(seurat_obj[["ATAC"]])))
cat(sprintf("Total dimensions: %d genes/peaks × %d cells\n", 
            nrow(seurat_obj[["RNA"]]) + nrow(seurat_obj[["ATAC"]]), 
            ncol(seurat_obj)))

# ============================================================================
# ATAC QUALITY CONTROL METRICS
# ============================================================================
cat("\n=== Computing ATAC QC metrics ===\n\n")

DefaultAssay(seurat_obj) <- "ATAC"

# Nucleosome signal
cat("Computing nucleosome signal...\n")
seurat_obj <- NucleosomeSignal(object = seurat_obj)
seurat_obj$nucleosome_group <- ifelse(
  seurat_obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4'
)

# TSS Enrichment
cat("Computing TSS enrichment scores...\n")
seurat_obj <- TSSEnrichment(seurat_obj)

# ============================================================================
# RNA QUALITY CONTROL METRICS
# ============================================================================
cat("\n=== Computing RNA QC metrics ===\n\n")

DefaultAssay(seurat_obj) <- "RNA"

# Calculate mitochondrial percentage
cat(sprintf("Computing mitochondrial percentage (pattern: %s)...\n", MT_PATTERN))
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = MT_PATTERN)

# ============================================================================
# PRINT PRE-FILTER SUMMARY
# ============================================================================
cat("\n=== Pre-filtering QC Summary ===\n\n")

qc_df <- seurat_obj@meta.data

cat("RNA metrics:\n")
cat(sprintf("  nCount_RNA:     median=%.0f, range=[%d, %d]\n",
            median(qc_df$nCount_RNA), min(qc_df$nCount_RNA), max(qc_df$nCount_RNA)))
cat(sprintf("  nFeature_RNA:   median=%.0f, range=[%d, %d]\n",
            median(qc_df$nFeature_RNA), min(qc_df$nFeature_RNA), max(qc_df$nFeature_RNA)))
cat(sprintf("  percent.mt:     median=%.1f%%, range=[%.1f%%, %.1f%%]\n",
            median(qc_df$percent.mt), min(qc_df$percent.mt), max(qc_df$percent.mt)))

cat("\nATAC metrics:\n")
cat(sprintf("  nCount_ATAC:        median=%.0f, range=[%d, %d]\n",
            median(qc_df$nCount_ATAC), min(qc_df$nCount_ATAC), max(qc_df$nCount_ATAC)))
cat(sprintf("  TSS.enrichment:     median=%.2f, range=[%.2f, %.2f]\n",
            median(qc_df$TSS.enrichment), min(qc_df$TSS.enrichment), max(qc_df$TSS.enrichment)))
cat(sprintf("  nucleosome_signal:  median=%.2f, range=[%.2f, %.2f]\n",
            median(qc_df$nucleosome_signal), min(qc_df$nucleosome_signal), max(qc_df$nucleosome_signal)))

# ============================================================================
# SAVE PRE-FILTER QC DATA FOR PLOTTING
# ============================================================================
# Store pre-filter data for before/after comparison
qc_df_before <- qc_df
qc_df_before$filter_status <- "Before Filtering"
n_cells_before <- ncol(seurat_obj)

# ============================================================================
# APPLY QC FILTERS
# ============================================================================
cat("\n=== Applying QC filters ===\n\n")

cat("Filter thresholds:\n")
cat(sprintf("  RNA features:        %d - %d\n", QC_MIN_FEATURES_RNA, QC_MAX_FEATURES_RNA))
cat(sprintf("  Mitochondrial %%:     < %d%%\n", QC_MAX_PERCENT_MT))
cat(sprintf("  ATAC counts:         %d - %d\n", QC_MIN_COUNT_ATAC, QC_MAX_COUNT_ATAC))
cat(sprintf("  RNA counts:          %d - %d\n", QC_MIN_COUNT_RNA, QC_MAX_COUNT_RNA))
cat(sprintf("  Nucleosome signal:   < %.1f\n", QC_MAX_NUCLEOSOME_SIGNAL))
cat(sprintf("  TSS enrichment:      > %d\n", QC_MIN_TSS_ENRICHMENT))

# Apply filters
seurat_obj <- subset(
  x = seurat_obj,
  subset = 
    nFeature_RNA > QC_MIN_FEATURES_RNA & 
    nFeature_RNA < QC_MAX_FEATURES_RNA & 
    percent.mt < QC_MAX_PERCENT_MT &
    nCount_ATAC < QC_MAX_COUNT_ATAC &
    nCount_ATAC > QC_MIN_COUNT_ATAC &
    nCount_RNA < QC_MAX_COUNT_RNA &
    nCount_RNA > QC_MIN_COUNT_RNA &
    nucleosome_signal < QC_MAX_NUCLEOSOME_SIGNAL &
    TSS.enrichment > QC_MIN_TSS_ENRICHMENT
)

n_cells_after <- ncol(seurat_obj)
cat(sprintf("\nCells before filtering: %d\n", n_cells_before))
cat(sprintf("Cells after filtering:  %d\n", n_cells_after))
cat(sprintf("Cells removed:          %d (%.1f%%)\n", 
            n_cells_before - n_cells_after,
            100 * (n_cells_before - n_cells_after) / n_cells_before))

# ============================================================================
# GENERATE QC VIOLIN PLOTS (BEFORE AND AFTER FILTERING)
# ============================================================================
cat("\n=== Generating QC violin plots ===\n\n")

# Create QC plots directory
qc_plots_dir <- file.path(OUTPUT_SEURAT_DIR, "qc_plots")
if (!dir.exists(qc_plots_dir)) {
  dir.create(qc_plots_dir, recursive = TRUE)
}

# ============================================================================
# FACETED VIOLIN PLOT (LIKE SEURAT VlnPlot) - PRE-FILTER
# ============================================================================
# Create a nice faceted violin plot similar to your manual plot

# Reshape data to long format for faceting
qc_metrics <- c("TSS.enrichment", "nCount_ATAC", "nCount_RNA", 
                "nFeature_RNA", "nucleosome_signal", "percent.mt")

qc_long_before <- qc_df_before %>%
  dplyr::select(all_of(qc_metrics)) %>%
  mutate(group = SAMPLE_NAME) %>%
  tidyr::pivot_longer(cols = all_of(qc_metrics), 
                      names_to = "metric", 
                      values_to = "value")

# Set factor levels for nice ordering
qc_long_before$metric <- factor(qc_long_before$metric, levels = qc_metrics)

# Create faceted violin plot - PRE-FILTER
facet_violin_before <- ggplot(qc_long_before, aes(x = group, y = value)) +
  geom_violin(fill = "white", color = "black", trim = TRUE, scale = "width") +
  geom_jitter(width = 0.3, size = 0.1, alpha = 0.3) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  labs(x = "group", y = "value", 
       title = sprintf("%s: QC Metrics (Before Filtering)", SAMPLE_NAME),
       subtitle = sprintf("Total cells: %d", n_cells_before)) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save pre-filter faceted violin
facet_before_file <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_facet_violin_BEFORE.pdf"))
ggsave(facet_before_file, facet_violin_before, width = 12, height = 8)
cat(sprintf("Saved faceted violin (before): %s\n", facet_before_file))

facet_before_png <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_facet_violin_BEFORE.png"))
ggsave(facet_before_png, facet_violin_before, width = 12, height = 8, dpi = 150)

# Get post-filter QC data
qc_df_after <- seurat_obj@meta.data
qc_df_after$filter_status <- "After Filtering"

# Create faceted violin plot - POST-FILTER
qc_long_after <- qc_df_after %>%
  dplyr::select(all_of(qc_metrics)) %>%
  mutate(group = SAMPLE_NAME) %>%
  tidyr::pivot_longer(cols = all_of(qc_metrics), 
                      names_to = "metric", 
                      values_to = "value")
qc_long_after$metric <- factor(qc_long_after$metric, levels = qc_metrics)

facet_violin_after <- ggplot(qc_long_after, aes(x = group, y = value)) +
  geom_violin(fill = "white", color = "black", trim = TRUE, scale = "width") +
  geom_jitter(width = 0.3, size = 0.1, alpha = 0.3) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  labs(x = "group", y = "value",
       title = sprintf("%s: QC Metrics (After Filtering)", SAMPLE_NAME),
       subtitle = sprintf("Filtered cells: %d (removed %d, %.1f%%)", 
                          n_cells_after, n_cells_before - n_cells_after,
                          100 * (n_cells_before - n_cells_after) / n_cells_before)) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save post-filter faceted violin
facet_after_file <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_facet_violin_AFTER.pdf"))
ggsave(facet_after_file, facet_violin_after, width = 12, height = 8)
cat(sprintf("Saved faceted violin (after): %s\n", facet_after_file))

facet_after_png <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_facet_violin_AFTER.png"))
ggsave(facet_after_png, facet_violin_after, width = 12, height = 8, dpi = 150)

# ============================================================================
# SIDE-BY-SIDE BEFORE/AFTER FACETED VIOLIN (COMBINED)
# ============================================================================
# Combine for side-by-side comparison
qc_long_before$status <- "Before"
qc_long_after$status <- "After"
qc_long_combined <- rbind(qc_long_before, qc_long_after)
qc_long_combined$status <- factor(qc_long_combined$status, levels = c("Before", "After"))

facet_violin_combined <- ggplot(qc_long_combined, aes(x = status, y = value, fill = status)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.05, alpha = 0.2) +
  scale_fill_manual(values = c("Before" = "#E74C3C", "After" = "#27AE60")) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(x = "", y = "value",
       title = sprintf("%s: QC Metrics Before vs After Filtering", SAMPLE_NAME),
       subtitle = sprintf("Before: %d cells | After: %d cells | Removed: %d (%.1f%%)",
                          n_cells_before, n_cells_after,
                          n_cells_before - n_cells_after,
                          100 * (n_cells_before - n_cells_after) / n_cells_before)) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

facet_combined_file <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_facet_violin_COMBINED.pdf"))
ggsave(facet_combined_file, facet_violin_combined, width = 12, height = 10)
cat(sprintf("Saved faceted violin (combined): %s\n", facet_combined_file))

facet_combined_png <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_facet_violin_COMBINED.png"))
ggsave(facet_combined_png, facet_violin_combined, width = 12, height = 10, dpi = 150)

# ============================================================================
# DETAILED COMPARISON VIOLIN PLOTS
# ============================================================================

# Combine before and after data for comparison plots
qc_combined <- rbind(
  qc_df_before[, c("nCount_RNA", "nFeature_RNA", "percent.mt", 
                   "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "filter_status")],
  qc_df_after[, c("nCount_RNA", "nFeature_RNA", "percent.mt", 
                  "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "filter_status")]
)
qc_combined$filter_status <- factor(qc_combined$filter_status, 
                                    levels = c("Before Filtering", "After Filtering"))

# Helper function to create violin plot with threshold lines
create_qc_violin <- function(data, y_var, y_label, threshold_low = NULL, threshold_high = NULL,
                             log_scale = FALSE) {
  p <- ggplot(data, aes(x = filter_status, y = .data[[y_var]], fill = filter_status)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.8) +
    scale_fill_manual(values = c("Before Filtering" = "#E74C3C", "After Filtering" = "#27AE60")) +
    labs(y = y_label, x = "", title = y_label) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Add threshold lines
  if (!is.null(threshold_low)) {
    p <- p + geom_hline(yintercept = threshold_low, linetype = "dashed", color = "blue", linewidth = 0.8)
  }
  if (!is.null(threshold_high)) {
    p <- p + geom_hline(yintercept = threshold_high, linetype = "dashed", color = "red", linewidth = 0.8)
  }
  
  if (log_scale) {
    p <- p + scale_y_log10()
  }
  
  return(p)
}

# Create individual plots with threshold lines
p_ncount_rna <- create_qc_violin(qc_combined, "nCount_RNA", "RNA Counts (nCount_RNA)",
                                  threshold_low = QC_MIN_COUNT_RNA, 
                                  threshold_high = QC_MAX_COUNT_RNA,
                                  log_scale = TRUE)

p_nfeature_rna <- create_qc_violin(qc_combined, "nFeature_RNA", "RNA Features (nFeature_RNA)",
                                    threshold_low = QC_MIN_FEATURES_RNA,
                                    threshold_high = QC_MAX_FEATURES_RNA)

p_percent_mt <- create_qc_violin(qc_combined, "percent.mt", "Mitochondrial %",
                                  threshold_high = QC_MAX_PERCENT_MT)

p_ncount_atac <- create_qc_violin(qc_combined, "nCount_ATAC", "ATAC Counts (nCount_ATAC)",
                                   threshold_low = QC_MIN_COUNT_ATAC,
                                   threshold_high = QC_MAX_COUNT_ATAC,
                                   log_scale = TRUE)

p_tss <- create_qc_violin(qc_combined, "TSS.enrichment", "TSS Enrichment",
                           threshold_low = QC_MIN_TSS_ENRICHMENT)

p_nucleosome <- create_qc_violin(qc_combined, "nucleosome_signal", "Nucleosome Signal",
                                  threshold_high = QC_MAX_NUCLEOSOME_SIGNAL)

# Combine all plots into one figure
combined_plot <- (p_ncount_rna | p_nfeature_rna | p_percent_mt) /
                 (p_ncount_atac | p_tss | p_nucleosome) +
  plot_annotation(
    title = sprintf("%s: QC Metrics Before and After Filtering", SAMPLE_NAME),
    subtitle = sprintf("Before: %d cells | After: %d cells | Removed: %d (%.1f%%)",
                       n_cells_before, n_cells_after, 
                       n_cells_before - n_cells_after,
                       100 * (n_cells_before - n_cells_after) / n_cells_before),
    caption = "Blue dashed = lower threshold | Red dashed = upper threshold",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 11))
  )

# Save combined plot
qc_plot_file <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_violin_before_after.pdf"))
ggsave(qc_plot_file, combined_plot, width = 14, height = 10)
cat(sprintf("Saved QC violin plots: %s\n", qc_plot_file))

# Also save as PNG for quick viewing
qc_plot_png <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_violin_before_after.png"))
ggsave(qc_plot_png, combined_plot, width = 14, height = 10, dpi = 150)
cat(sprintf("Saved QC violin plots: %s\n", qc_plot_png))

# ============================================================================
# GENERATE DISTRIBUTION HISTOGRAMS (PRE-FILTER ONLY FOR THRESHOLD SELECTION)
# ============================================================================
cat("\n=== Generating pre-filter distribution histograms ===\n\n")

# Create histograms with current threshold markers (helps decide new thresholds)
create_histogram <- function(data, x_var, x_label, threshold_low = NULL, threshold_high = NULL,
                             bins = 50, log_scale = FALSE) {
  p <- ggplot(data, aes(x = .data[[x_var]])) +
    geom_histogram(bins = bins, fill = "#3498DB", color = "black", alpha = 0.7) +
    labs(x = x_label, y = "Cell Count", title = paste("Distribution:", x_label)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
  
  if (!is.null(threshold_low)) {
    p <- p + geom_vline(xintercept = threshold_low, linetype = "dashed", color = "blue", linewidth = 1) +
      annotate("text", x = threshold_low, y = Inf, label = sprintf("Min: %s", threshold_low),
               vjust = 2, hjust = -0.1, color = "blue", size = 3)
  }
  if (!is.null(threshold_high)) {
    p <- p + geom_vline(xintercept = threshold_high, linetype = "dashed", color = "red", linewidth = 1) +
      annotate("text", x = threshold_high, y = Inf, label = sprintf("Max: %s", threshold_high),
               vjust = 2, hjust = 1.1, color = "red", size = 3)
  }
  
  if (log_scale) {
    p <- p + scale_x_log10()
  }
  
  return(p)
}

# Use pre-filter data for histogram
h_ncount_rna <- create_histogram(qc_df_before, "nCount_RNA", "RNA Counts",
                                  threshold_low = QC_MIN_COUNT_RNA, 
                                  threshold_high = QC_MAX_COUNT_RNA,
                                  log_scale = TRUE)

h_nfeature_rna <- create_histogram(qc_df_before, "nFeature_RNA", "RNA Features",
                                    threshold_low = QC_MIN_FEATURES_RNA,
                                    threshold_high = QC_MAX_FEATURES_RNA)

h_percent_mt <- create_histogram(qc_df_before, "percent.mt", "Mitochondrial %",
                                  threshold_high = QC_MAX_PERCENT_MT)

h_ncount_atac <- create_histogram(qc_df_before, "nCount_ATAC", "ATAC Counts",
                                   threshold_low = QC_MIN_COUNT_ATAC,
                                   threshold_high = QC_MAX_COUNT_ATAC,
                                   log_scale = TRUE)

h_tss <- create_histogram(qc_df_before, "TSS.enrichment", "TSS Enrichment",
                           threshold_low = QC_MIN_TSS_ENRICHMENT)

h_nucleosome <- create_histogram(qc_df_before, "nucleosome_signal", "Nucleosome Signal",
                                  threshold_high = QC_MAX_NUCLEOSOME_SIGNAL)

# Combine histograms
hist_plot <- (h_ncount_rna | h_nfeature_rna | h_percent_mt) /
             (h_ncount_atac | h_tss | h_nucleosome) +
  plot_annotation(
    title = sprintf("%s: Pre-Filter QC Distributions with Current Thresholds", SAMPLE_NAME),
    subtitle = sprintf("Total cells: %d | Use these plots to adjust QC thresholds in config.R", n_cells_before),
    caption = "Blue = lower threshold | Red = upper threshold",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 11))
  )

# Save histogram plots
hist_plot_file <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_histograms_prefilter.pdf"))
ggsave(hist_plot_file, hist_plot, width = 14, height = 10)
cat(sprintf("Saved pre-filter histograms: %s\n", hist_plot_file))

hist_plot_png <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_histograms_prefilter.png"))
ggsave(hist_plot_png, hist_plot, width = 14, height = 10, dpi = 150)
cat(sprintf("Saved pre-filter histograms: %s\n", hist_plot_png))

# ============================================================================
# SCATTER PLOTS FOR QC METRICS
# ============================================================================
cat("\n=== Generating QC scatter plots ===\n\n")

# RNA scatter: nCount vs nFeature colored by MT%
scatter_rna <- ggplot(qc_df_before, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_color_viridis_c(option = "plasma", name = "MT %") +
  geom_hline(yintercept = c(QC_MIN_FEATURES_RNA, QC_MAX_FEATURES_RNA), 
             linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(QC_MIN_COUNT_RNA, QC_MAX_COUNT_RNA), 
             linetype = "dashed", color = "blue") +
  scale_x_log10() +
  labs(x = "RNA Counts (log10)", y = "RNA Features", 
       title = "RNA: Counts vs Features (colored by MT%)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11))

# ATAC scatter: nCount vs TSS colored by nucleosome signal
scatter_atac <- ggplot(qc_df_before, aes(x = nCount_ATAC, y = TSS.enrichment, color = nucleosome_signal)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_color_viridis_c(option = "viridis", name = "Nucleosome\nSignal") +
  geom_hline(yintercept = QC_MIN_TSS_ENRICHMENT, linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(QC_MIN_COUNT_ATAC, QC_MAX_COUNT_ATAC), 
             linetype = "dashed", color = "blue") +
  scale_x_log10() +
  labs(x = "ATAC Counts (log10)", y = "TSS Enrichment",
       title = "ATAC: Counts vs TSS (colored by Nucleosome Signal)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11))

# MT% vs Features
scatter_mt <- ggplot(qc_df_before, aes(x = percent.mt, y = nFeature_RNA)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#3498DB") +
  geom_vline(xintercept = QC_MAX_PERCENT_MT, linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(QC_MIN_FEATURES_RNA, QC_MAX_FEATURES_RNA), 
             linetype = "dashed", color = "blue") +
  labs(x = "Mitochondrial %", y = "RNA Features",
       title = "Mitochondrial % vs RNA Features") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11))

# Nucleosome vs TSS
scatter_nuc <- ggplot(qc_df_before, aes(x = nucleosome_signal, y = TSS.enrichment)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#E74C3C") +
  geom_vline(xintercept = QC_MAX_NUCLEOSOME_SIGNAL, linetype = "dashed", color = "red") +
  geom_hline(yintercept = QC_MIN_TSS_ENRICHMENT, linetype = "dashed", color = "blue") +
  labs(x = "Nucleosome Signal", y = "TSS Enrichment",
       title = "Nucleosome Signal vs TSS Enrichment") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11))

# Combine scatter plots
scatter_plot <- (scatter_rna | scatter_atac) / (scatter_mt | scatter_nuc) +
  plot_annotation(
    title = sprintf("%s: QC Scatter Plots (Pre-Filter)", SAMPLE_NAME),
    subtitle = "Dashed lines show current filter thresholds",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 11))
  )

scatter_plot_file <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_scatter_prefilter.pdf"))
ggsave(scatter_plot_file, scatter_plot, width = 14, height = 10)
cat(sprintf("Saved scatter plots: %s\n", scatter_plot_file))

scatter_plot_png <- file.path(qc_plots_dir, paste0(SAMPLE_NAME, "_qc_scatter_prefilter.png"))
ggsave(scatter_plot_png, scatter_plot, width = 14, height = 10, dpi = 150)
cat(sprintf("Saved scatter plots: %s\n", scatter_plot_png))

# ============================================================================
# FILTER LOW-EXPRESSION GENES
# ============================================================================
cat("\n=== Filtering low-expression genes ===\n\n")

DefaultAssay(seurat_obj) <- "RNA"

# Number of cells in the dataset
num_cells <- ncol(seurat_obj)

# Calculate the fraction of cells each gene is expressed in
gene_counts <- GetAssayData(seurat_obj, layer = "counts")
gene_frac <- Matrix::rowSums(gene_counts > 0) / num_cells

# Keep genes expressed in at least the specified fraction of cells
n_genes_before <- nrow(seurat_obj[["RNA"]])
genes_keep <- names(gene_frac[gene_frac >= QC_MIN_GENE_FRACTION])

# Subset only the RNA assay to preserve ATAC assay
seurat_obj[["RNA"]] <- subset(seurat_obj[["RNA"]], features = genes_keep)

n_genes_after <- nrow(seurat_obj[["RNA"]])

cat(sprintf("Gene filter: Keep genes expressed in >= %.0f%% of cells\n", 
            QC_MIN_GENE_FRACTION * 100))
cat(sprintf("Genes before: %d\n", n_genes_before))
cat(sprintf("Genes after:  %d\n", n_genes_after))
cat(sprintf("Genes removed: %d (%.1f%%)\n",
            n_genes_before - n_genes_after,
            100 * (n_genes_before - n_genes_after) / n_genes_before))

# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n=== Final dataset summary ===\n\n")

cat(sprintf("Assays: %s\n", paste(Assays(seurat_obj), collapse = ", ")))
cat(sprintf("RNA genes:    %d\n", nrow(seurat_obj[["RNA"]])))
cat(sprintf("ATAC peaks:   %d\n", nrow(seurat_obj[["ATAC"]])))
cat(sprintf("Total cells:  %d\n", ncol(seurat_obj)))

# ============================================================================
# SAVE OUTPUT
# ============================================================================
cat("\n=== Saving output ===\n\n")

output_file <- file.path(
  OUTPUT_SEURAT_DIR, 
  paste0(SAMPLE_NAME, "_seurat_multiome_processed.rds")
)

saveRDS(seurat_obj, output_file)
cat(sprintf("Saved: %s\n", output_file))

# Save QC summary as CSV for reference
qc_summary <- data.frame(
  sample = SAMPLE_NAME,
  cells_raw = n_cells_before,
  cells_filtered = n_cells_after,
  cells_removed = n_cells_before - n_cells_after,
  genes_raw = n_genes_before,
  genes_filtered = n_genes_after,
  peaks = nrow(seurat_obj[["ATAC"]]),
  qc_min_features_rna = QC_MIN_FEATURES_RNA,
  qc_max_features_rna = QC_MAX_FEATURES_RNA,
  qc_max_percent_mt = QC_MAX_PERCENT_MT,
  qc_min_count_atac = QC_MIN_COUNT_ATAC,
  qc_max_count_atac = QC_MAX_COUNT_ATAC,
  qc_min_tss = QC_MIN_TSS_ENRICHMENT,
  qc_max_nucleosome = QC_MAX_NUCLEOSOME_SIGNAL,
  processing_date = Sys.time()
)

qc_file <- file.path(OUTPUT_SEURAT_DIR, paste0(SAMPLE_NAME, "_qc_summary.csv"))
write.csv(qc_summary, qc_file, row.names = FALSE)
cat(sprintf("Saved QC summary: %s\n", qc_file))

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 031 COMPLETE\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Output files:\n")
cat(sprintf("  - %s\n", basename(output_file)))
cat(sprintf("  - %s\n", basename(qc_file)))
cat("\nQC plots (in qc_plots/ subdirectory):\n")
cat(sprintf("  - %s_qc_facet_violin_BEFORE.pdf/png   (Faceted violin - pre-filter)\n", SAMPLE_NAME))
cat(sprintf("  - %s_qc_facet_violin_AFTER.pdf/png    (Faceted violin - post-filter)\n", SAMPLE_NAME))
cat(sprintf("  - %s_qc_facet_violin_COMBINED.pdf/png (Before/After side-by-side)\n", SAMPLE_NAME))
cat(sprintf("  - %s_qc_violin_before_after.pdf/png   (Detailed with thresholds)\n", SAMPLE_NAME))
cat(sprintf("  - %s_qc_histograms_prefilter.pdf/png  (Distributions)\n", SAMPLE_NAME))
cat(sprintf("  - %s_qc_scatter_prefilter.pdf/png     (2D scatter plots)\n", SAMPLE_NAME))
cat("\n*** IMPORTANT: Review the QC plots to adjust thresholds in config.R if needed ***\n")

cat("\nNext step: Run 02a_data_splitting.R\n")
