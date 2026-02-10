# ============================================================================
# Step_040: Gene-Specific Feature Extraction (Automated)
# ============================================================================
# This script extracts ATAC peaks within gene windows and adds HVGs as features for prediction.
# It reads configuration from config.R.
#
# APPROACH:
# - For each target gene (HVG), define a genomic window (e.g., ±250kb from TSS)
# - Extract all ATAC peaks that fall within this window
# - Create feature matrices (cells × peaks × HVGs) for train/validation/test splits
# - Save organized feature data for efficient model training
#
# Input: Smoothed data from Step_030
# Output: Gene-specific feature matrices for model training
#
# Usage:
#   1. Ensure config.R is properly configured
#   2. Run: Rscript 04_feature_extraction.R
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
  library(Signac)
  library(GenomicRanges)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(reshape2)
  library(viridis)
})

# Load species-specific annotation packages
load_annotation_packages()

# ============================================================================
# PLOTTING THEME
# ============================================================================
theme_pub <- function(base_size = 12, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Title and subtitle
      plot.title = element_text(size = base_size * 1.2, face = "bold", 
                                hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size * 0.9, hjust = 0.5, 
                                   color = "grey30", margin = margin(b = 10)),
      plot.caption = element_text(size = base_size * 0.7, hjust = 1, 
                                  color = "grey50"),
      
      # Axis
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size * 0.9, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      
      # Panel
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white"),
      
      # Legend
      legend.title = element_text(size = base_size * 0.9, face = "bold"),
      legend.text = element_text(size = base_size * 0.8),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.position = "right",
      
      # Facets
      strip.text = element_text(size = base_size * 0.9, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      
      # Margins
      plot.margin = margin(15, 15, 15, 15)
    )
}

# Colorblind-friendly palette (Wong 2011)
colorblind_palette <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7",  # reddish purple
  "#999999"   # grey
)

# Colors for data splits
split_colors_pub <- c(
  "train" = "#0072B2",       # blue
  "validation" = "#E69F00", # orange
  "test" = "#009E73"        # green
)

# Helper function to save plots in both PNG and PDF
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 600) {
  # Save PNG (raster, high resolution)
  ggsave(paste0(filename, ".png"), plot, width = width, height = height, dpi = dpi)
  
  # Save PDF (vector, publication quality)
  cairo_pdf(paste0(filename, ".pdf"), width = width, height = height)
  print(plot)
  dev.off()
}

# Create output directories
create_output_directories()

# Set parameters from config
gene_window_kb <- GENE_WINDOW_KB
min_peaks_per_gene <- MIN_PEAKS_PER_GENE
n_hvg_genes <- N_HVG_GENES
n_hvg_features <- N_HVG_FEATURES
seed <- SEED_FEATURES

set.seed(seed)

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 060: Feature Extraction\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Parameters:\n")
cat(sprintf("  Gene window: ±%d kb from TSS\n", gene_window_kb))
cat(sprintf("  Min peaks per gene: %d\n", min_peaks_per_gene))
cat(sprintf("  HVG target genes: %d\n", n_hvg_genes))
cat(sprintf("  HVG features: %d\n", n_hvg_features))
cat(sprintf("  Random seed: %d\n\n", seed))

# ============================================================================
# LOAD PREPROCESSED DATA
# ============================================================================
cat("=== Loading preprocessed data ===\n\n")

# Load original Seurat object (for annotations)
seurat_file <- file.path(
  OUTPUT_SPLITS_DIR, 
  paste0(SAMPLE_NAME, "_seurat_obj_with_splits.rds")
)

if (!file.exists(seurat_file)) {
  stop(sprintf("ERROR: Seurat file not found: %s\nPlease run Step_040 first.", seurat_file))
}

cat(sprintf("Loading Seurat object: %s\n", seurat_file))
seurat_obj <- readRDS(seurat_file)
cat(sprintf("  Loaded: %d cells, %d genes, %d peaks\n\n",
            ncol(seurat_obj), nrow(seurat_obj[["RNA"]]), nrow(seurat_obj[["ATAC"]])))

# Load smoothed data for all splits
cat("Loading smoothed data...\n")
smoothed_train <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_train.rds"))
smoothed_val <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_validation.rds"))
smoothed_test <- readRDS(file.path(OUTPUT_METACELLS_DIR, "smoothed_test.rds"))

cat(sprintf("  Train: %d cells\n", length(smoothed_train$cells)))
cat(sprintf("  Validation: %d cells\n", length(smoothed_val$cells)))
cat(sprintf("  Test: %d cells\n", length(smoothed_test$cells)))

# ============================================================================
# GET GENE ANNOTATIONS
# ============================================================================
cat("\n=== Extracting gene annotations ===\n\n")

ensdb <- get_ensdb()
gene_annotation <- GetGRangesFromEnsDb(ensdb = ensdb)

# Convert to UCSC chromosome naming
if (SPECIES == "mouse") {
  seqlevels(gene_annotation) <- paste0('chr', seqlevels(gene_annotation))
} else {
  seqlevelsStyle(gene_annotation) <- "UCSC"
}

cat(sprintf("Gene annotations: %d entries\n", length(gene_annotation)))

# ============================================================================
# FIND HIGHLY VARIABLE GENES (FROM TRAINING ONLY)
# ============================================================================
cat("\n=== Finding highly variable genes from TRAINING data only ===\n\n")

train_cells_for_hvg <- colnames(seurat_obj)[seurat_obj$data_split == "train"]
seurat_train_subset <- seurat_obj[, train_cells_for_hvg]

cat(sprintf("Using %d training cells (excluding %d val + %d test cells)\n",
            length(train_cells_for_hvg),
            sum(seurat_obj$data_split == "validation"),
            sum(seurat_obj$data_split == "test")))

DefaultAssay(seurat_train_subset) <- "RNA"

cat("Computing highly variable genes from training data...\n")
seurat_train_subset <- FindVariableFeatures(
  seurat_train_subset,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

all_hvgs <- VariableFeatures(seurat_train_subset)
cat(sprintf("✓ Computed %d highly variable genes (HVGs) from training data only\n", length(all_hvgs)))
cat("✓ No data leakage: validation/test data NOT used in HVG selection\n")

# ============================================================================
# FIND VARIABLE ATAC PEAKS (FROM TRAINING ONLY)
# ============================================================================
cat("\n=== Finding variable ATAC peaks from TRAINING data only ===\n\n")

DefaultAssay(seurat_train_subset) <- "ATAC"

cat("Computing top variable ATAC peaks from training data...\n")
seurat_train_subset <- FindTopFeatures(
  seurat_train_subset,
  assay = "ATAC",
  min.cutoff = "q5",
  verbose = FALSE
)

variable_peaks <- VariableFeatures(seurat_train_subset)
cat(sprintf("✓ Computed %d variable ATAC peaks from training data only\n", length(variable_peaks)))

# Fallback if no variable peaks found
if (length(variable_peaks) == 0) {
  cat("WARNING: FindTopFeatures returned 0 peaks. Trying alternative approach...\n")
  seurat_train_subset <- FindTopFeatures(seurat_train_subset, assay = "ATAC", min.cutoff = "q0", verbose = FALSE)
  variable_peaks <- VariableFeatures(seurat_train_subset)
  cat(sprintf("After retry: %d variable ATAC peaks\n", length(variable_peaks)))
}

cat("✓ No data leakage: validation/test data NOT used in peak selection\n")

# ============================================================================
# DEFINE TOP HVG FEATURES
# ============================================================================
cat("\n=== Defining top HVG features ===\n\n")

top_hvg_features <- head(all_hvgs, min(n_hvg_features, length(all_hvgs)))
cat(sprintf("Selected top %d HVGs as additional expression features\n", length(top_hvg_features)))
cat("First 10 HVG features:", paste(head(top_hvg_features, 10), collapse = ", "), "\n")

# ============================================================================
# DEFINE TARGET GENE SETS
# ============================================================================
cat("\n=== Defining target gene sets ===\n\n")

# ============================================================================
# DEFINE TARGET GENE SETS
# ============================================================================
cat("\n=== Defining target gene sets ===\n\n")

# Gene Set 1: Use all genes from HVG_GENE_FILE only
hvg_genes <- NULL
if (!is.null(HVG_GENE_FILE) && HVG_GENE_FILE != "" && file.exists(HVG_GENE_FILE)) {
  hvg_genes <- readLines(HVG_GENE_FILE)
  hvg_genes <- hvg_genes[nchar(hvg_genes) > 0]
  cat(sprintf("Loaded %d HVG genes from user-supplied file: %s\n", length(hvg_genes), HVG_GENE_FILE))
  cat("First 10 HVG genes:", paste(head(hvg_genes, 10), collapse = ", "), "\n")
} else if (!is.null(HVG_GENE_FILE) && HVG_GENE_FILE != "") {
  stop(sprintf("HVG_GENE_FILE is set but file does not exist: %s", HVG_GENE_FILE))
}

# Gene Set 2: Use all genes from RANDOM_GENE_FILE only (do not fall back)
random_genes <- NULL
if (!is.null(RANDOM_GENE_FILE) && RANDOM_GENE_FILE != "" && file.exists(RANDOM_GENE_FILE)) {
  random_genes <- readLines(RANDOM_GENE_FILE)
  random_genes <- random_genes[nchar(random_genes) > 0]
  cat(sprintf("Loaded %d random genes from user-supplied file: %s\n", length(random_genes), RANDOM_GENE_FILE))
  cat("First 10 random genes:", paste(head(random_genes, 10), collapse = ", "), "\n")
} else if (!is.null(RANDOM_GENE_FILE) && RANDOM_GENE_FILE != "") {
  stop(sprintf("RANDOM_GENE_FILE is set but file does not exist: %s", RANDOM_GENE_FILE))
}

# Create list of gene sets to process
gene_sets <- list()
if (!is.null(hvg_genes) && length(hvg_genes) > 0) {
  gene_sets$HVG <- hvg_genes
}
if (!is.null(random_genes) && length(random_genes) > 0) {
  gene_sets$Random_genes <- random_genes
}

cat(sprintf("\nTotal gene sets to process: %d\n", length(gene_sets)))
for (set_name in names(gene_sets)) {
  cat(sprintf("  - %s: %d genes\n", set_name, length(gene_sets[[set_name]])))
}

# ============================================================================
# GET PEAK ANNOTATIONS
# ============================================================================
cat("\n=== Extracting peak annotations ===\n\n")

DefaultAssay(seurat_obj) <- "ATAC"
peak_ranges <- granges(seurat_obj[["ATAC"]])

cat(sprintf("Total ATAC peaks: %d\n", length(peak_ranges)))
cat("\nFirst 5 peaks:\n")
print(head(peak_ranges, 5))

# ============================================================================
# HELPER FUNCTIONS FOR FEATURE EXTRACTION
# ============================================================================

# Extract peak features for a specific gene
extract_gene_peak_features <- function(gene_name, peak_indices, smoothed_data) {
  all_peak_names <- rownames(smoothed_data$atac_log1p)
  gene_peaks <- all_peak_names[peak_indices]
  feature_matrix <- smoothed_data$atac_log1p[gene_peaks, , drop = FALSE]
  t(feature_matrix)  # Return cells × peaks
}

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
  t(hvg_matrix)  # Return cells × genes
}

# Extract target gene expression
extract_gene_expression <- function(gene_name, smoothed_data) {
  all_gene_names <- rownames(smoothed_data$rna_log1p)
  
  if (!gene_name %in% all_gene_names) {
    warning(sprintf("Gene %s not found in RNA data", gene_name))
    return(NULL)
  }
  as.numeric(smoothed_data$rna_log1p[gene_name, ])
}

# ============================================================================
# PROCESS EACH GENE SET
# ============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("=== PROCESSING GENE SETS ===\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

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
  
  # Filter gene annotations to target genes
  cat(sprintf("\nFiltering gene annotations to target genes...\n"))
  target_gene_idx <- which(gene_annotation$gene_name %in% target_genes)
  gene_annotation_set <- gene_annotation[target_gene_idx]
  
  cat(sprintf("  Target genes in this set: %d\n", length(target_genes)))
  cat(sprintf("  Matched: %d annotation entries\n", length(gene_annotation_set)))
  
  # Collapse multiple transcripts per gene
  cat("  Collapsing multiple transcripts per gene...\n")
  unique_gene_names <- unique(gene_annotation_set$gene_name)
  gene_annotation_set <- gene_annotation_set[match(unique_gene_names, gene_annotation_set$gene_name)]
  
  cat(sprintf("  Target genes with annotations: %d\n", length(gene_annotation_set)))
  
  # Check for missing genes
  missing_genes <- setdiff(target_genes, gene_annotation_set$gene_name)
  if (length(missing_genes) > 0) {
    cat(sprintf("  Warning: %d target genes lack annotations:\n", length(missing_genes)))
    cat("    ", paste(head(missing_genes, 10), collapse = ", "), "\n")
    if (length(missing_genes) > 10) cat("    ... and", length(missing_genes) - 10, "more\n")
  }
  
  if (length(gene_annotation_set) == 0) {
    cat(sprintf("ERROR: No genes in %s set have annotations. Skipping.\n", set_name))
    next
  }
  
  # ===== Define gene windows and extract peaks =====
  cat(sprintf("\n=== Defining gene windows (±%dkb from TSS) ===\n", gene_window_kb))
  
  window_bp <- gene_window_kb * 1000
  gene_windows <- resize(gene_annotation_set, width = 1, fix = "start")
  gene_windows <- promoters(gene_windows, upstream = window_bp, downstream = window_bp)
  
  cat(sprintf("Created genomic windows for %d genes\n", length(gene_windows)))
  
  # Filter to variable peaks
  cat("\nFiltering peaks to variable features only...\n")
  peak_names <- names(peak_ranges)
  variable_peak_indices <- which(peak_names %in% variable_peaks)
  
  if (length(variable_peak_indices) == 0) {
    cat("  WARNING: No variable peaks found! Using all peaks instead.\n")
    peak_ranges_variable <- peak_ranges
    variable_peak_indices <- seq_along(peak_ranges)
  } else {
    peak_ranges_variable <- peak_ranges[variable_peak_indices]
  }
  
  cat(sprintf("  Total peaks: %d\n", length(peak_ranges)))
  cat(sprintf("  Variable peaks used: %d (%.1f%%)\n", 
              length(peak_ranges_variable), 
              100 * length(peak_ranges_variable) / length(peak_ranges)))
  
  # Find overlaps between peaks and gene windows
  cat("\nFinding variable peaks within gene windows...\n")
  overlaps <- findOverlaps(peak_ranges_variable, gene_windows)
  
  gene_to_peaks <- split(variable_peak_indices[queryHits(overlaps)], subjectHits(overlaps))
  names(gene_to_peaks) <- gene_windows$gene_name[as.integer(names(gene_to_peaks))]
  
  # Filter genes with too few peaks
  peaks_per_gene <- sapply(gene_to_peaks, length)
  
  cat(sprintf("\nPeaks per gene statistics:\n"))
  cat(sprintf("  Mean: %.1f\n", mean(peaks_per_gene)))
  cat(sprintf("  Median: %.0f\n", median(peaks_per_gene)))
  cat(sprintf("  Min: %d\n", min(peaks_per_gene)))
  cat(sprintf("  Max: %d\n", max(peaks_per_gene)))
  
  valid_genes <- names(gene_to_peaks)[peaks_per_gene >= min_peaks_per_gene]
  gene_to_peaks <- gene_to_peaks[valid_genes]
  
  cat(sprintf("\nGenes with ≥%d peaks: %d (%.1f%% of annotated targets)\n",
              min_peaks_per_gene, length(gene_to_peaks),
              100 * length(gene_to_peaks) / length(gene_annotation_set)))
  
  # ===== Extract features for each gene =====
  cat("\n=== Extracting peak features for each gene ===\n")
  cat("This may take a few minutes...\n\n")
  
  gene_features <- list()
  n_genes <- length(gene_to_peaks)
  progress_interval <- max(1, floor(n_genes / 20))
  
  for (i in seq_along(gene_to_peaks)) {
    gene_name <- names(gene_to_peaks)[i]
    peak_indices <- gene_to_peaks[[i]]
    
    if (i %% progress_interval == 0 || i == n_genes) {
      cat(sprintf("  [%d/%d] %.1f%% complete\r", i, n_genes, 100 * i / n_genes))
    }
    
    tryCatch({
      # Extract gene-specific peak features
      train_peaks <- extract_gene_peak_features(gene_name, peak_indices, smoothed_train)
      val_peaks <- extract_gene_peak_features(gene_name, peak_indices, smoothed_val)
      test_peaks <- extract_gene_peak_features(gene_name, peak_indices, smoothed_test)
      
      # Extract HVG expression features
      train_hvgs <- extract_hvg_features(gene_name, top_hvg_features, smoothed_train)
      val_hvgs <- extract_hvg_features(gene_name, top_hvg_features, smoothed_val)
      test_hvgs <- extract_hvg_features(gene_name, top_hvg_features, smoothed_test)
      
      # Combine features
      train_X <- cbind(train_peaks, train_hvgs)
      val_X <- cbind(val_peaks, val_hvgs)
      test_X <- cbind(test_peaks, test_hvgs)
      
      # Extract target gene expression
      train_y <- extract_gene_expression(gene_name, smoothed_train)
      val_y <- extract_gene_expression(gene_name, smoothed_val)
      test_y <- extract_gene_expression(gene_name, smoothed_test)
      
      # Store in organized structure
      gene_features[[gene_name]] <- list(
        gene_name = gene_name,
        n_peaks = length(peak_indices),
        n_hvg_features = ncol(train_hvgs),
        n_total_features = ncol(train_X),
        peak_names = rownames(smoothed_train$atac_log1p)[peak_indices],
        hvg_feature_names = colnames(train_hvgs),
        
        train = list(X = train_X, y = train_y, cells = smoothed_train$cells),
        validation = list(X = val_X, y = val_y, cells = smoothed_val$cells),
        test = list(X = test_X, y = test_y, cells = smoothed_test$cells)
      )
    }, error = function(e) {
      cat(sprintf("\nWarning: Failed to extract features for gene %s: %s\n", gene_name, e$message))
    })
  }
  
  cat("\n")
  
  # ===== Quality control and summary =====
  cat("\n=== Feature Extraction Summary for", set_name, "===\n")
  
  gene_features <- gene_features[!sapply(gene_features, is.null)]
  
  cat(sprintf("Successfully extracted features for %d genes\n", length(gene_features)))
  
  if (length(gene_features) > 0) {
    n_peaks_per_gene <- sapply(gene_features, function(x) x$n_peaks)
    n_hvg_per_gene <- sapply(gene_features, function(x) x$n_hvg_features)
    n_total_per_gene <- sapply(gene_features, function(x) x$n_total_features)
    
    cat(sprintf("\nPeaks per gene:\n"))
    cat(sprintf("  Mean: %.1f, Median: %.0f, Range: [%d, %d]\n", 
                mean(n_peaks_per_gene), median(n_peaks_per_gene), 
                min(n_peaks_per_gene), max(n_peaks_per_gene)))
    
    cat(sprintf("\nTotal features per gene:\n"))
    cat(sprintf("  Mean: %.1f, Median: %.0f, Range: [%d, %d]\n", 
                mean(n_total_per_gene), median(n_total_per_gene), 
                min(n_total_per_gene), max(n_total_per_gene)))
    
    # Example output
    first_gene <- gene_features[[1]]
    cat(sprintf("\nExample: Gene '%s'\n", first_gene$gene_name))
    cat(sprintf("  Peak features: %d, HVG features: %d, Total: %d\n",
                first_gene$n_peaks, first_gene$n_hvg_features, first_gene$n_total_features))
    cat(sprintf("  Train: %d × %d, Val: %d × %d, Test: %d × %d\n",
                nrow(first_gene$train$X), ncol(first_gene$train$X),
                nrow(first_gene$validation$X), ncol(first_gene$validation$X),
                nrow(first_gene$test$X), ncol(first_gene$test$X)))
  }
  
  # ===== Save outputs =====
  cat("\n=== Saving extracted features ===\n")
  
  # Save complete feature list
  output_file <- file.path(set_output_dir, "gene_specific_features.rds")
  saveRDS(gene_features, output_file)
  cat(sprintf("Saved: %s\n", output_file))
  
  # Save metadata
  if (length(gene_features) > 0) {
    feature_metadata <- data.frame(
      gene_name = names(gene_features),
      n_peaks = sapply(gene_features, function(x) x$n_peaks),
      n_hvg_features = sapply(gene_features, function(x) x$n_hvg_features),
      n_total_features = sapply(gene_features, function(x) x$n_total_features),
      n_train_cells = sapply(gene_features, function(x) nrow(x$train$X)),
      n_val_cells = sapply(gene_features, function(x) nrow(x$validation$X)),
      n_test_cells = sapply(gene_features, function(x) nrow(x$test$X)),
      stringsAsFactors = FALSE
    )
    
    metadata_file <- file.path(set_output_dir, "gene_features_metadata.csv")
    write.csv(feature_metadata, metadata_file, row.names = FALSE)
    cat(sprintf("Saved metadata: %s\n", metadata_file))
  }
  
  # Save parameters
  params <- list(
    gene_set = set_name,
    sample_name = SAMPLE_NAME,
    gene_window_kb = gene_window_kb,
    min_peaks_per_gene = min_peaks_per_gene,
    n_hvg_features = n_hvg_features,
    n_genes_extracted = length(gene_features),
    seed = seed,
    extraction_date = Sys.time()
  )
  
  params_file <- file.path(set_output_dir, "feature_extraction_params.rds")
  saveRDS(params, params_file)
  cat(sprintf("Saved parameters: %s\n", params_file))
  
  cat(sprintf("\n>>> Completed %s gene set <<<\n", set_name))
  cat(paste(rep("-", 70), collapse = ""), "\n")
}

# ============================================================================
# DIAGNOSTIC PLOTS
# ============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("=== GENERATING DIAGNOSTIC PLOTS ===\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Create plots directory
plots_dir <- file.path(OUTPUT_FEATURES_DIR, "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}
cat("Plots directory:", plots_dir, "\n\n")

# Load all gene sets for comparison
all_gene_features <- list()
for (set_name in names(gene_sets)) {
  set_dir <- file.path(OUTPUT_FEATURES_DIR, set_name)
  feature_file <- file.path(set_dir, "gene_specific_features.rds")
  if (file.exists(feature_file)) {
    all_gene_features[[set_name]] <- readRDS(feature_file)
    cat(sprintf("Loaded %s gene set: %d genes\n", set_name, length(all_gene_features[[set_name]])))
  }
}

# Use primary gene set for detailed plots
main_set <- names(gene_sets)[1]
gene_features <- all_gene_features[[main_set]]

cat("\nGenerating plots based on all", length(all_gene_features), "gene sets\n")
cat("Primary gene set for detailed plots:", main_set, "(", length(gene_features), "genes)\n\n")

# ----------------------------------------------------------------------------
# PLOT 1: Distribution of peaks per target gene (ALL GENE SETS)
# ----------------------------------------------------------------------------
cat("Creating Plot 1: Peaks per target gene distribution (all gene sets)...\n")

# Combine data from all gene sets
peaks_all_sets <- lapply(names(all_gene_features), function(set_name) {
  gf <- all_gene_features[[set_name]]
  data.frame(
    gene = names(gf),
    n_peaks = sapply(gf, function(x) x$n_peaks),
    gene_set = set_name
  )
})
peaks_all_df <- do.call(rbind, peaks_all_sets)

# Color palette for gene sets
set_colors <- colorblind_palette[1:length(all_gene_features)]
names(set_colors) <- names(all_gene_features)

p1 <- ggplot(peaks_all_df, aes(x = n_peaks, fill = gene_set)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, 
                 alpha = 0.6, color = "white", position = "identity") +
  geom_density(aes(color = gene_set), linewidth = 1.2, fill = NA) +
  scale_fill_manual(values = set_colors, name = "Gene Set") +
  scale_color_manual(values = set_colors, name = "Gene Set") +
  labs(
    title = sprintf("Distribution of ATAC Peaks per Target Gene (±%d kb)", gene_window_kb),
    subtitle = sprintf("Comparing %d gene sets | Min peaks/gene ≥ %d", 
                       length(all_gene_features), min_peaks_per_gene),
    x = "Number of ATAC peaks within gene window",
    y = "Density"
  ) +
  theme_pub() +
  theme(legend.position = c(0.85, 0.85)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

save_plot(p1, file.path(plots_dir, "01_peaks_per_gene_distribution"), 
                      width = 9, height = 6)
cat("  ✓ Saved peaks per gene distribution (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 2: Total features distribution (ALL GENE SETS) - NEW PLOT
# ----------------------------------------------------------------------------
cat("Creating Plot 2: Total features per gene distribution (all gene sets)...\n")

# Combine total features from all gene sets
features_all_sets <- lapply(names(all_gene_features), function(set_name) {
  gf <- all_gene_features[[set_name]]
  data.frame(
    gene = names(gf),
    n_total = sapply(gf, function(x) x$n_total_features),
    n_peaks = sapply(gf, function(x) x$n_peaks),
    n_hvg = sapply(gf, function(x) x$n_hvg_features),
    gene_set = set_name
  )
})
features_all_df <- do.call(rbind, features_all_sets)

p2 <- ggplot(features_all_df, aes(x = n_total, fill = gene_set)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, 
                 alpha = 0.6, color = "white", position = "identity") +
  geom_density(aes(color = gene_set), linewidth = 1.2, fill = NA) +
  scale_fill_manual(values = set_colors, name = "Gene Set") +
  scale_color_manual(values = set_colors, name = "Gene Set") +
  geom_vline(data = aggregate(n_total ~ gene_set, features_all_df, median),
             aes(xintercept = n_total, color = gene_set), 
             linetype = "dashed", linewidth = 1) +
  labs(
    title = "Distribution of Total Features per Target Gene",
    subtitle = sprintf("Peak features + %d HVG expression features | All gene sets", n_hvg_features),
    x = "Total number of features (peaks + HVG)",
    y = "Density"
  ) +
  theme_pub() +
  theme(legend.position = c(0.85, 0.85)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

save_plot(p2, file.path(plots_dir, "02_total_features_distribution_comparison"), 
                      width = 9, height = 6)
cat("  ✓ Saved total features comparison (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 2a-2n: Individual total features distribution for EACH gene set
# ----------------------------------------------------------------------------
cat("Creating individual total features plots for each gene set...\n")

individual_feature_plots <- list()
for (i in seq_along(names(all_gene_features))) {
  set_name <- names(all_gene_features)[i]
  set_data <- features_all_df[features_all_df$gene_set == set_name, ]
  
  plot_letter <- letters[i]  # a, b, c, ...
  
  p_individual <- ggplot(set_data, aes(x = n_total)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, 
                   fill = set_colors[set_name], color = "white", alpha = 0.8) +
    geom_density(color = "black", linewidth = 1.2) +
    geom_vline(xintercept = median(set_data$n_total), 
               linetype = "dashed", color = colorblind_palette[6], linewidth = 1) +
    geom_vline(xintercept = mean(set_data$n_total), 
               linetype = "solid", color = colorblind_palette[6], linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
             label = sprintf("n = %d genes\nMean: %.1f\nMedian: %.0f\nMin: %d\nMax: %d", 
                             nrow(set_data),
                             mean(set_data$n_total),
                             median(set_data$n_total),
                             min(set_data$n_total),
                             max(set_data$n_total)),
             size = 3.5, fontface = "bold", color = "grey30") +
    labs(
      title = sprintf("Total Features per Gene: %s", set_name),
      subtitle = sprintf("Peak features (var) + %d HVG expression features", n_hvg_features),
      x = "Total number of features (peaks + HVG)",
      y = "Density"
    ) +
    theme_pub() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  individual_feature_plots[[set_name]] <- p_individual
  
  save_plot(p_individual, 
            file.path(plots_dir, sprintf("02%s_total_features_%s", plot_letter, set_name)), 
            width = 8, height = 6)
  cat(sprintf("  ✓ Saved total features for %s (PNG + PDF)\n", set_name))
}

# ----------------------------------------------------------------------------
# PLOT 3: Feature breakdown boxplot (ALL GENE SETS)
# ----------------------------------------------------------------------------
cat("Creating Plot 3: Feature breakdown by gene set...\n")

feature_breakdown <- reshape2::melt(features_all_df, 
                                     id.vars = c("gene", "gene_set"),
                                     measure.vars = c("n_peaks", "n_hvg"),
                                     variable.name = "feature_type",
                                     value.name = "count")
feature_breakdown$feature_type <- factor(feature_breakdown$feature_type,
                                          levels = c("n_peaks", "n_hvg"),
                                          labels = c("ATAC Peaks", "HVG Features"))

p3 <- ggplot(feature_breakdown, aes(x = gene_set, y = count, fill = feature_type)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.3, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("ATAC Peaks" = colorblind_palette[5], 
                                "HVG Features" = colorblind_palette[3]),
                    name = "Feature Type") +
  labs(
    title = "Feature Breakdown by Gene Set",
    subtitle = "Comparing peak features and HVG expression features",
    x = "Gene Set",
    y = "Number of Features"
  ) +
  theme_pub() +
  theme(legend.position = "top")

save_plot(p3, file.path(plots_dir, "03_feature_breakdown_by_set"), 
          width = 8, height = 6)
cat("  ✓ Saved feature breakdown (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 4: Cumulative distribution of peaks per gene (ALL GENE SETS)
# ----------------------------------------------------------------------------
cat("Creating Plot 4: Cumulative peaks distribution...\n")

peaks_cumulative <- lapply(names(all_gene_features), function(set_name) {
  df <- peaks_all_df[peaks_all_df$gene_set == set_name, ]
  df <- df[order(df$n_peaks), ]
  df$cumulative <- seq_len(nrow(df)) / nrow(df)
  df
})
peaks_cumulative_df <- do.call(rbind, peaks_cumulative)

p4 <- ggplot(peaks_cumulative_df, aes(x = n_peaks, y = cumulative, color = gene_set)) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = set_colors, name = "Gene Set") +
  labs(
    title = "Cumulative Distribution of Peaks per Gene",
    subtitle = "All gene sets compared",
    x = "Number of ATAC peaks",
    y = "Cumulative proportion of genes"
  ) +
  theme_pub() +
  theme(legend.position = c(0.85, 0.25)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1))

save_plot(p4, file.path(plots_dir, "04_peaks_cumulative_distribution"), 
          width = 8, height = 6)
cat("  ✓ Saved cumulative distribution (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 5: Distance from peaks to gene TSS (PRIMARY SET)
# ----------------------------------------------------------------------------
cat("Creating Plot 5: Peak distance to TSS...\n")

# Calculate distances for a sample of genes
set.seed(seed)
sample_genes <- sample(names(gene_features), min(100, length(gene_features)))

distance_data <- lapply(sample_genes, function(gene_name) {
  gf <- gene_features[[gene_name]]
  peak_names <- gf$peak_names
  
  # Get gene TSS
  gene_idx <- which(gene_annotation$gene_name == gene_name)
  if (length(gene_idx) == 0) return(NULL)
  
  gene_gr <- gene_annotation[gene_idx[1]]
  tss <- ifelse(as.character(strand(gene_gr)) == "+", start(gene_gr), end(gene_gr))
  gene_chr <- as.character(seqnames(gene_gr))
  
  # Parse peak positions
  distances <- sapply(peak_names, function(pn) {
    parts <- strsplit(pn, "-")[[1]]
    if (length(parts) >= 3) {
      peak_chr <- parts[1]
      peak_start <- as.numeric(parts[2])
      peak_end <- as.numeric(parts[3])
      peak_mid <- (peak_start + peak_end) / 2
      
      if (peak_chr == gene_chr) {
        return(peak_mid - tss)
      }
    }
    return(NA)
  })
  
  data.frame(gene = gene_name, distance_kb = distances / 1000)
})

distance_df <- do.call(rbind, distance_data)
distance_df <- distance_df[!is.na(distance_df$distance_kb), ]

p5 <- ggplot(distance_df, aes(x = distance_kb)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100, 
                 fill = colorblind_palette[2], color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "solid", color = colorblind_palette[6], linewidth = 1) +
  annotate("text", x = 0, y = Inf, vjust = 2, hjust = -0.1,
           label = "TSS", color = colorblind_palette[6], fontface = "bold", size = 4) +
  labs(
    title = "Distance from ATAC Peaks to Gene TSS",
    subtitle = sprintf("Window: ±%d kb | Sample of %d genes, %d peaks", 
                       gene_window_kb, length(sample_genes), nrow(distance_df)),
    x = "Distance to TSS (kb)",
    y = "Density"
  ) +
  theme_pub() +
  scale_x_continuous(limits = c(-gene_window_kb - 10, gene_window_kb + 10))

save_plot(p5, file.path(plots_dir, "05_peak_distance_to_TSS"), 
          width = 8, height = 6)
cat("  ✓ Saved peak distance to TSS (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 6: Gene expression distribution by split
# ----------------------------------------------------------------------------
cat("Creating Plot 6: Target gene expression distribution...\n")

# Sample gene expression values
expr_data <- lapply(sample_genes, function(gene_name) {
  gf <- gene_features[[gene_name]]
  data.frame(
    gene = gene_name,
    expression = c(gf$train$y, gf$validation$y, gf$test$y),
    split = c(rep("train", length(gf$train$y)),
              rep("validation", length(gf$validation$y)),
              rep("test", length(gf$test$y)))
  )
})
expr_df <- do.call(rbind, expr_data)

p6 <- ggplot(expr_df, aes(x = expression, fill = split)) +
  geom_density(alpha = 0.6, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = split_colors_pub, name = "Data Split") +
  labs(
    title = "Distribution of Target Gene Expression by Data Split",
    subtitle = sprintf("Smoothed log1p expression | %d sampled genes", length(sample_genes)),
    x = "log1p(Expression)",
    y = "Density"
  ) +
  theme_pub() +
  theme(legend.position = c(0.85, 0.85))

save_plot(p6, file.path(plots_dir, "06_expression_distribution_by_split"), 
          width = 8, height = 6)
cat("  ✓ Saved expression distribution by split (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 7: Peak accessibility distribution
# ----------------------------------------------------------------------------
cat("Creating Plot 7: Peak accessibility distribution...\n")

# Get peak accessibility for sample genes
peak_data <- lapply(sample_genes[1:min(20, length(sample_genes))], function(gene_name) {
  gf <- gene_features[[gene_name]]
  peak_vals <- as.vector(gf$train$X[, 1:min(gf$n_peaks, 50)])
  data.frame(gene = gene_name, accessibility = peak_vals)
})
peak_df <- do.call(rbind, peak_data)

p7 <- ggplot(peak_df, aes(x = accessibility)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, 
                 fill = colorblind_palette[3], color = "white", alpha = 0.8) +
  geom_density(color = colorblind_palette[6], linewidth = 1.2) +
  labs(
    title = "Distribution of ATAC Peak Accessibility (Features)",
    subtitle = "Smoothed log1p accessibility from training data",
    x = "log1p(Accessibility)",
    y = "Density"
  ) +
  theme_pub()

save_plot(p7, file.path(plots_dir, "07_peak_accessibility_distribution"), 
          width = 8, height = 6)
cat("  ✓ Saved peak accessibility distribution (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 8: Feature composition per gene
# ----------------------------------------------------------------------------
cat("Creating Plot 8: Feature composition per gene...\n")

feature_comp <- data.frame(
  gene = names(gene_features),
  peaks = sapply(gene_features, function(x) x$n_peaks),
  hvg = sapply(gene_features, function(x) x$n_hvg_features)
)

feature_comp_long <- reshape2::melt(feature_comp, id.vars = "gene", 
                                     variable.name = "feature_type", 
                                     value.name = "count")
feature_comp_long$feature_type <- factor(feature_comp_long$feature_type, 
                                          levels = c("peaks", "hvg"),
                                          labels = c("ATAC Peaks", "HVG Expression"))

# Order genes by total features
gene_order <- feature_comp$gene[order(feature_comp$peaks + feature_comp$hvg)]
feature_comp_long$gene <- factor(feature_comp_long$gene, levels = gene_order)

# Show only top and bottom genes
n_show <- min(30, length(gene_features))
genes_to_show <- c(head(gene_order, n_show/2), tail(gene_order, n_show/2))
feature_comp_subset <- feature_comp_long[feature_comp_long$gene %in% genes_to_show, ]

p8 <- ggplot(feature_comp_subset, aes(x = gene, y = count, fill = feature_type)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("ATAC Peaks" = colorblind_palette[5], 
                                "HVG Expression" = colorblind_palette[1]),
                    name = "Feature Type") +
  labs(
    title = "Feature Composition per Target Gene",
    subtitle = sprintf("Showing %d genes with lowest and highest feature counts", n_show),
    x = "Target Gene",
    y = "Number of Features"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "top")

save_plot(p8, file.path(plots_dir, "08_feature_composition_per_gene"), 
          width = 10, height = 6)
cat("  ✓ Saved feature composition (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 9: Peak-Expression correlation for sample genes
# ----------------------------------------------------------------------------
cat("Creating Plot 9: Peak-expression correlation...\n")

# Calculate correlation between peaks and expression for sample genes
cor_data <- lapply(sample_genes[1:min(50, length(sample_genes))], function(gene_name) {
  gf <- gene_features[[gene_name]]
  n_peaks_use <- min(gf$n_peaks, 20)
  
  cors <- sapply(1:n_peaks_use, function(i) {
    cor(gf$train$X[, i], gf$train$y, use = "complete.obs")
  })
  
  data.frame(
    gene = gene_name,
    peak_idx = 1:n_peaks_use,
    correlation = cors
  )
})
cor_df <- do.call(rbind, cor_data)
cor_df <- cor_df[!is.na(cor_df$correlation), ]

p9 <- ggplot(cor_df, aes(x = correlation)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, 
                 fill = colorblind_palette[7], color = "white", alpha = 0.8) +
  geom_density(color = "black", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = mean(cor_df$correlation), 
             linetype = "solid", color = colorblind_palette[6], linewidth = 1) +
  annotate("text", x = mean(cor_df$correlation), y = Inf, vjust = 2, hjust = -0.1,
           label = sprintf("Mean: %.3f", mean(cor_df$correlation)),
           color = colorblind_palette[6], fontface = "bold", size = 4) +
  labs(
    title = "Correlation Between Peak Accessibility and Gene Expression",
    subtitle = sprintf("Training data | %d genes, %d peak-gene pairs", 
                       length(unique(cor_df$gene)), nrow(cor_df)),
    x = "Pearson Correlation",
    y = "Density"
  ) +
  theme_pub() +
  scale_x_continuous(limits = c(-0.5, 0.5))

save_plot(p9, file.path(plots_dir, "09_peak_expression_correlation"), 
          width = 8, height = 6)
cat("  ✓ Saved peak-expression correlation (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 10: Cells per split bar chart
# ----------------------------------------------------------------------------
cat("Creating Plot 10: Cells per split...\n")

split_counts <- data.frame(
  split = c("train", "validation", "test"),
  cells = c(length(smoothed_train$cells), 
            length(smoothed_val$cells), 
            length(smoothed_test$cells))
)
split_counts$split <- factor(split_counts$split, levels = c("train", "validation", "test"))

p10 <- ggplot(split_counts, aes(x = split, y = cells, fill = split)) +
  geom_col(alpha = 0.9, color = "black", linewidth = 0.5) +
  geom_text(aes(label = format(cells, big.mark = ",")), 
            vjust = -0.5, fontface = "bold", size = 4) +
  scale_fill_manual(values = split_colors_pub) +
  labs(
    title = "Number of Cells per Data Split",
    subtitle = sprintf("Total: %s cells", 
                       format(sum(split_counts$cells), big.mark = ",")),
    x = "Data Split",
    y = "Number of Cells"
  ) +
  theme_pub() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

save_plot(p10, file.path(plots_dir, "10_cells_per_split"), 
          width = 6, height = 6)
cat("  ✓ Saved cells per split (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 11: Feature dimensions boxplot (all gene sets)
# ----------------------------------------------------------------------------
cat("Creating Plot 11: Feature dimensions summary (all gene sets)...\n")

feature_dims_all <- lapply(names(all_gene_features), function(set_name) {
  gf <- all_gene_features[[set_name]]
  data.frame(
    gene_set = set_name,
    type = rep(c("Peak Features", "HVG Features", "Total Features"), 
               each = length(gf)),
    count = c(
      sapply(gf, function(x) x$n_peaks),
      sapply(gf, function(x) x$n_hvg_features),
      sapply(gf, function(x) x$n_total_features)
    )
  )
})
feature_dims_df <- do.call(rbind, feature_dims_all)
feature_dims_df$type <- factor(feature_dims_df$type, 
                                levels = c("Peak Features", "HVG Features", "Total Features"))

p11 <- ggplot(feature_dims_df, aes(x = type, y = count, fill = gene_set)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.3, position = position_dodge(0.8)) +
  scale_fill_manual(values = set_colors, name = "Gene Set") +
  labs(
    title = "Distribution of Feature Dimensions by Gene Set",
    subtitle = sprintf("Comparing %d gene sets", length(all_gene_features)),
    x = "Feature Type",
    y = "Number of Features"
  ) +
  theme_pub() +
  theme(legend.position = "top")

save_plot(p11, file.path(plots_dir, "11_feature_dimensions_by_set"), 
          width = 9, height = 6)
cat("  ✓ Saved feature dimensions boxplot (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 12: Top peaks per gene heatmap (for a few genes)
# ----------------------------------------------------------------------------
cat("Creating Plot 12: Peak accessibility heatmap...\n")

# Select top 8 genes with most peaks
top_genes <- names(sort(sapply(gene_features, function(x) x$n_peaks), decreasing = TRUE)[1:min(8, length(gene_features))])

# Create heatmap data
heatmap_data <- lapply(top_genes, function(gene_name) {
  gf <- gene_features[[gene_name]]
  n_peaks_use <- min(gf$n_peaks, 15)  # Limit to 15 peaks per gene
  
  mean_access <- colMeans(gf$train$X[, 1:n_peaks_use, drop = FALSE])
  data.frame(
    gene = gene_name,
    peak = paste0("Peak_", 1:n_peaks_use),
    accessibility = mean_access
  )
})
heatmap_df <- do.call(rbind, heatmap_data)

p12 <- ggplot(heatmap_df, aes(x = peak, y = gene, fill = accessibility)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis(option = "plasma", name = "Mean\nAccessibility") +
  labs(
    title = "Mean Peak Accessibility for Top Target Genes",
    subtitle = "Training data | Top 8 genes by peak count (up to 15 peaks each)",
    x = "Peak Index (ordered by position)",
    y = "Target Gene"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9))

save_plot(p12, file.path(plots_dir, "12_peak_accessibility_heatmap"), 
          width = 10, height = 6)
cat("  ✓ Saved peak accessibility heatmap (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 13: Gene set comparison summary
# ----------------------------------------------------------------------------
cat("Creating Plot 13: Gene set comparison...\n")

set_summary <- data.frame(
  gene_set = names(all_gene_features),
  n_genes = sapply(all_gene_features, length),
  median_peaks = sapply(all_gene_features, function(gf) median(sapply(gf, function(x) x$n_peaks))),
  median_total = sapply(all_gene_features, function(gf) median(sapply(gf, function(x) x$n_total_features)))
)

set_summary_long <- reshape2::melt(set_summary, id.vars = "gene_set",
                                    variable.name = "metric", 
                                    value.name = "value")
set_summary_long$metric <- factor(set_summary_long$metric,
                                   levels = c("n_genes", "median_peaks", "median_total"),
                                   labels = c("Number of Genes", "Median Peaks/Gene", "Median Total Features"))

p13 <- ggplot(set_summary_long, aes(x = gene_set, y = value, fill = gene_set)) +
  geom_col(alpha = 0.9, color = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f", value)), 
            vjust = -0.3, fontface = "bold", size = 3.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = set_colors) +
  labs(
    title = "Gene Set Comparison Summary",
    subtitle = "Key metrics for each target gene set",
    x = "Gene Set",
    y = "Value"
  ) +
  theme_pub() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

save_plot(p13, file.path(plots_dir, "13_gene_set_comparison"), 
          width = 10, height = 5)
cat("  ✓ Saved gene set comparison (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# PLOT 14: Summary statistics panel
# ----------------------------------------------------------------------------
cat("Creating Plot 14: Summary panel...\n")

# Build summary for all gene sets
set_summary_text <- paste(sapply(names(all_gene_features), function(set_name) {
  gf <- all_gene_features[[set_name]]
  sprintf("  %s: %d genes, median %.1f peaks/gene", 
          set_name, length(gf), median(sapply(gf, function(x) x$n_peaks)))
}), collapse = "\n")

summary_text <- paste0(
  "FEATURE EXTRACTION SUMMARY\n\n",
  sprintf("Sample: %s\n\n", SAMPLE_NAME),
  "Parameters:\n",
  sprintf("  Gene window: ±%d kb from TSS\n", gene_window_kb),
  sprintf("  Min peaks/gene: %d\n", min_peaks_per_gene),
  sprintf("  HVG features: %d\n\n", n_hvg_features),
  "Gene Sets:\n",
  set_summary_text, "\n\n",
  "Data Dimensions:\n",
  sprintf("  Train cells: %d\n", length(smoothed_train$cells)),
  sprintf("  Validation cells: %d\n", length(smoothed_val$cells)),
  sprintf("  Test cells: %d\n\n", length(smoothed_test$cells)),
  sprintf("Primary Set (%s) Statistics:\n", main_set),
  sprintf("  Median peaks/gene: %.1f\n", median(sapply(gene_features, function(x) x$n_peaks))),
  sprintf("  Median total features/gene: %.1f\n", median(sapply(gene_features, function(x) x$n_total_features)))
)

p14 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = summary_text, 
           hjust = 0.5, vjust = 0.5, size = 3.5, family = "mono") +
  theme_void() +
  theme(plot.margin = margin(20, 20, 20, 20),
        plot.background = element_rect(fill = "white", color = NA))

save_plot(p14, file.path(plots_dir, "14_summary_panel"), 
                      width = 7, height = 7)
cat("  ✓ Saved summary panel (PNG + PDF)\n")

# ----------------------------------------------------------------------------
# COMBINED PDF REPORT
# ----------------------------------------------------------------------------
cat("\nCreating combined PDF report...\n")

cairo_pdf(file.path(plots_dir, "Step_040_Feature_Extraction_Report.pdf"), 
          width = 10, height = 8)

# Title page
plot.new()
text(0.5, 0.6, "Step 040: Feature Extraction Report", cex = 2, font = 2)
text(0.5, 0.45, sprintf("Sample: %s", SAMPLE_NAME), cex = 1.4)
text(0.5, 0.35, sprintf("Date: %s", Sys.Date()), cex = 1.2)
text(0.5, 0.25, "Figures", cex = 1, font = 3)
text(0.5, 0.15, sprintf("Gene Sets: %s", paste(names(all_gene_features), collapse = ", ")), cex = 1)

# Print all plots
print(p1)
print(p2)
# Print individual feature plots for each gene set
for (set_name in names(individual_feature_plots)) {
  print(individual_feature_plots[[set_name]])
}
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
print(p10)
print(p11)
print(p12)
print(p13)
print(p14)

dev.off()
cat("  ✓ Saved combined PDF report (vector graphics)\n")

cat("\n")
cat(" plots saved to:", plots_dir, "\n")
cat("  - Each plot saved as both .png (600 dpi) and .pdf (vector)\n")
cat("  - Combined report: Step_040_Feature_Extraction_Report.pdf\n")

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 040 COMPLETE\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Gene-specific ATAC peak features extracted successfully\n\n")

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
cat(sprintf("  - Sample: %s\n", SAMPLE_NAME))
cat(sprintf("  - Gene window: ±%d kb from TSS\n", gene_window_kb))
cat(sprintf("  - Min peaks per gene: %d\n", min_peaks_per_gene))
cat(sprintf("  - HVG features: %d\n", n_hvg_features))
cat(sprintf("  - Train cells: %d\n", length(smoothed_train$cells)))
cat(sprintf("  - Validation cells: %d\n", length(smoothed_val$cells)))
cat(sprintf("  - Test cells: %d\n", length(smoothed_test$cells)))

cat("\nData leakage prevention:\n")
cat("  ✓ HVGs computed from TRAINING data only\n")
cat("  ✓ Variable ATAC peaks computed from TRAINING data only\n")
cat("  ✓ Target gene excluded from HVG features\n")
cat("  ✓ Feature selection does NOT use validation/test information\n")

cat("\nPublication-ready plots (in plots/ subdirectory, PNG + PDF formats):\n")
cat("  01 - Peaks per gene distribution (all gene sets comparison)\n")
cat("  02 - Total features per gene distribution (all gene sets comparison)\n")
# List individual gene set plots
for (i in seq_along(names(all_gene_features))) {
  set_name <- names(all_gene_features)[i]
  cat(sprintf("  02%s - Total features distribution: %s (individual)\n", letters[i], set_name))
}
cat("  03 - Feature breakdown by gene set (boxplot)\n")
cat("  04 - Cumulative distribution of peaks per gene\n")
cat("  05 - Distance from ATAC peaks to gene TSS\n")
cat("  06 - Target gene expression distribution by split\n")
cat("  07 - Peak accessibility distribution\n")
cat("  08 - Feature composition per gene (peaks vs HVG)\n")
cat("  09 - Peak-expression correlation distribution\n")
cat("  10 - Cells per data split\n")
cat("  11 - Feature dimensions by gene set (boxplot)\n")
cat("  12 - Peak accessibility heatmap (top genes)\n")
cat("  13 - Gene set comparison summary\n")
cat("  14 - Summary statistics panel\n")
cat("  + Combined PDF report (vector graphics)\n")

cat("\nOutput directory:", OUTPUT_FEATURES_DIR, "\n\n")

cat("Next step: Run model training scripts\n")
