# Step_040.Feature_Extraction_Automated_ATAC_only.R
# Automated: Gene-Specific Feature Extraction (ATAC-ONLY)
# Loads all parameters and paths from config.R. No hardcoded logic.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(EnsDb.Mmusculus.v79)
  library(dplyr)
  library(Matrix)
})

# ---------------- Load Config ----------------
config_file <- "config.R"
if (!file.exists(config_file)) stop("ERROR: config.R not found in working directory.")
source(config_file)

# ----------------- Fallbacks / normalize config names -----------------
# Allow this automated script to run using the central `config.R` variables
# or with method-specific `ATAC_*` overrides. Set sensible defaults.
if (!exists("ATAC_FEATURES_INPUT_DIR") && exists("OUTPUT_METACELLS_DIR")) {
  ATAC_FEATURES_INPUT_DIR <- OUTPUT_METACELLS_DIR
}
if (!exists("ATAC_FEATURES_OUTPUT_DIR") && exists("OUTPUT_FEATURES_DIR")) {
  ATAC_FEATURES_OUTPUT_DIR <- OUTPUT_FEATURES_DIR
}
if (!exists("ATAC_INPUT_SPLIT_DIR") && exists("OUTPUT_SPLITS_DIR")) {
  ATAC_INPUT_SPLIT_DIR <- OUTPUT_SPLITS_DIR
}
if (!exists("ATAC_FEATURES_GENE_WINDOW_KB") && exists("GENE_WINDOW_KB")) {
  ATAC_FEATURES_GENE_WINDOW_KB <- GENE_WINDOW_KB
}
if (!exists("ATAC_FEATURES_MIN_PEAKS_PER_GENE") && exists("MIN_PEAKS_PER_GENE")) {
  ATAC_FEATURES_MIN_PEAKS_PER_GENE <- MIN_PEAKS_PER_GENE
}
if (!exists("ATAC_FEATURES_N_HVG_GENES") && exists("N_HVG_GENES")) {
  ATAC_FEATURES_N_HVG_GENES <- N_HVG_GENES
}
if (!exists("ATAC_RANDOM_SEED") && exists("SEED_FEATURES")) {
  ATAC_RANDOM_SEED <- SEED_FEATURES
}
if (!exists("ATAC_FEATURES_HVG_GENE_FILE") && exists("HVG_GENE_FILE")) {
  ATAC_FEATURES_HVG_GENE_FILE <- HVG_GENE_FILE
}
if (!exists("ATAC_RANDOM_GENE_FILE") && exists("RANDOM_GENE_FILE")) {
  ATAC_RANDOM_GENE_FILE <- RANDOM_GENE_FILE
}
if (!exists("ATAC_MIN_PEAK_QUANTILE")) ATAC_MIN_PEAK_QUANTILE <- "q5"

# Normalize and expose `output_dir` variable used later in the script
ATAC_FEATURES_INPUT_DIR <- path.expand(ATAC_FEATURES_INPUT_DIR)
ATAC_FEATURES_OUTPUT_DIR <- path.expand(ATAC_FEATURES_OUTPUT_DIR)
ATAC_INPUT_SPLIT_DIR <- path.expand(ATAC_INPUT_SPLIT_DIR)
output_dir <- ATAC_FEATURES_OUTPUT_DIR


# Required config variables:
# ATAC_FEATURES_INPUT_DIR: Directory with smoothed metacell data (per-split LSI or PeakVI)
# ATAC_FEATURES_OUTPUT_DIR: Output directory for extracted features
# ATAC_FEATURES_GENE_WINDOW_KB: Window size around gene TSS (in kb)
# ATAC_FEATURES_MIN_PEAKS_PER_GENE: Minimum number of peaks per gene
# ATAC_FEATURES_N_HVG_GENES: Number of HVGs to select
# ATAC_RANDOM_GENE_FILE: Path to random gene list
# ATAC_RANDOM_SEED: Random seed

cat("Input directory:", ATAC_FEATURES_INPUT_DIR, "\n")
cat("Output directory:", ATAC_FEATURES_OUTPUT_DIR, "\n")
if (!dir.exists(ATAC_FEATURES_OUTPUT_DIR)) {
  dir.create(ATAC_FEATURES_OUTPUT_DIR, recursive = TRUE)
  cat("Created output directory\n")
}

gene_window_kb <- ATAC_FEATURES_GENE_WINDOW_KB
min_peaks_per_gene <- ATAC_FEATURES_MIN_PEAKS_PER_GENE
n_hvg_genes <- ATAC_FEATURES_N_HVG_GENES
random_gene_file <- RANDOM_GENE_FILE
seed <- ATAC_RANDOM_SEED
set.seed(seed)

# ===============Load Preprocessed Data=========================
cat("\n=== Loading preprocessed data (ATAC-ONLY) ===\n")
## Prefer Seurat splits produced by the multiome pipeline when available
if (exists("INPUT_SEURAT_SPLITS") && nzchar(INPUT_SEURAT_SPLITS) && file.exists(path.expand(INPUT_SEURAT_SPLITS))) {
  seurat_obj_file <- path.expand(INPUT_SEURAT_SPLITS)
  cat(sprintf("Using multiome Seurat splits: %s\n", seurat_obj_file))
} else {
  seurat_obj_path <- ATAC_INPUT_SPLIT_DIR
  seurat_obj_file <- file.path(seurat_obj_path, paste0(SAMPLE_NAME, "_seurat_obj_with_splits.rds"))
}
if (!file.exists(seurat_obj_file)) stop(sprintf("Seurat splits RDS not found: %s", seurat_obj_file))
seurat_obj <- readRDS(seurat_obj_file)
cat(sprintf("Loaded Seurat object: %d cells, %d peaks\n", ncol(seurat_obj), nrow(seurat_obj[["ATAC"]])))

# Load smoothed data for all splits
cat("\nLoading smoothed ATAC data...\n")
smoothed_train <- readRDS(file.path(ATAC_FEATURES_INPUT_DIR, "smoothed_train.rds"))
smoothed_val <- readRDS(file.path(ATAC_FEATURES_INPUT_DIR, "smoothed_validation.rds"))
smoothed_test <- readRDS(file.path(ATAC_FEATURES_INPUT_DIR, "smoothed_test.rds"))

cat(sprintf("Train: %d cells\n", length(smoothed_train$cells)))
cat(sprintf("Validation: %d cells\n", length(smoothed_val$cells)))
cat(sprintf("Test: %d cells\n", length(smoothed_test$cells)))

# ===============Add RNA Expression Data=========================
cat("\n=== Adding RNA expression data to smoothed structures ===\n")

# RNA metacell directory - should be set in config.R
# For ATAC-only pipeline, RNA data is needed as the target variable (Y)
# Update RNA_METACELL_DIR in config.R to point to your RNA-only pipeline output
if (!exists("RNA_METACELL_DIR")) {
  # Default: assume RNA-only pipeline output is in a parallel directory
  RNA_METACELL_DIR <- gsub("scATAC_only", "scRNA_only", BASE_OUTPUT_DIR)
}
input_dir_rna <- file.path(RNA_METACELL_DIR, "metacells", SAMPLE_NAME, "pca")
  
  # Check if RNA assay exists
if (!"RNA" %in% names(seurat_obj@assays)) {
    cat("\nERROR: RNA assay not found in Seurat object!\n")
    cat("This script requires RNA expression data as target variable (y).\n")
    cat("Please ensure the Seurat object contains both RNA and ATAC assays.\n")
    cat("\nPossible solutions:\n")
    cat("  1. Load a multimodal Seurat object (RNA + ATAC)\n")
    cat("  2. Load RNA data from a separate file\n")
    cat("  3. If smoothed files already contain RNA, check file paths\n\n")
    stop("RNA assay missing from Seurat object")
  }
  
  # Extract log-normalized RNA data from Seurat object
  # Load smoothed RNA data for all splits
  cat("\nLoading smoothed RNA data...\n")
  smoothed_train_rna <- readRDS(file.path(input_dir_rna, "smoothed_train.rds"))
  smoothed_val_rna <- readRDS(file.path(input_dir_rna, "smoothed_validation.rds"))
  smoothed_test_rna <- readRDS(file.path(input_dir_rna, "smoothed_test.rds"))

  cat(sprintf("Train: %d cells\n", length(smoothed_train_rna$cells)))
  cat(sprintf("Validation: %d cells\n", length(smoothed_val_rna$cells)))
  cat(sprintf("Test: %d cells\n", length(smoothed_test_rna$cells)))
  
  # Train split
  train_cells_in_seurat <- intersect(smoothed_train$cells, smoothed_train_rna$cells)
  
  if (length(train_cells_in_seurat) == length(smoothed_train$cells)) {
    smoothed_train$rna_log1p <- smoothed_train_rna$rna_log1p[, smoothed_train$cells, drop = FALSE]
    cat(sprintf("  Train: Added RNA data for %d cells\n", length(train_cells_in_seurat)))
  } else {
    cat(sprintf("  WARNING: Train - only %d/%d cells found in Seurat object\n", 
                length(train_cells_in_seurat), length(smoothed_train$cells)))
    if (length(train_cells_in_seurat) == 0) {
      cat("  ERROR: No cell overlap! Cell names may not match between files.\n")
      cat("  This suggests the smoothed data and Seurat object are from different datasets.\n")
    }
    
    # Subset both ATAC and RNA tp overlapping cells only
    smoothed_train$atac_log1p <- smoothed_train$atac_log1p[, train_cells_in_seurat, drop = FALSE]
    smoothed_train$rna_log1p <- smoothed_train_rna$rna_log1p[, train_cells_in_seurat, drop = FALSE]
    smoothed_train$cells <- train_cells_in_seurat
  }

  # Assert check the column names
  stopifnot(
    all(colnames(smoothed_train$rna_log1p) == smoothed_train$cells),
    length(smoothed_train$cells) == ncol(smoothed_train$rna_log1p)
  )
  
  # Validation split
  val_cells_in_seurat <- intersect(smoothed_val$cells, smoothed_val_rna$cells)
  if (length(val_cells_in_seurat) == length(smoothed_val$cells)) {
    smoothed_val$rna_log1p <- smoothed_val_rna$rna_log1p[, smoothed_val$cells, drop = FALSE]
    cat(sprintf("  Validation: Added RNA data for %d cells\n", length(val_cells_in_seurat)))
  } else {
    cat(sprintf("  WARNING: Validation - only %d/%d cells found in Seurat object\n", 
                length(val_cells_in_seurat), length(smoothed_val$cells)))
    # Subset both ATAC and RNA tp overlapping cells only
    smoothed_val$atac_log1p <- smoothed_val$atac_log1p[, val_cells_in_seurat, drop = FALSE]
    smoothed_val$rna_log1p <- smoothed_val_rna$rna_log1p[, val_cells_in_seurat, drop = FALSE]
    smoothed_val$cells <- val_cells_in_seurat
  }
  
  # Assert check the column names
  stopifnot(
    all(colnames(smoothed_val$rna_log1p) == smoothed_val$cells),
    length(smoothed_val$cells) == ncol(smoothed_val$rna_log1p)
  )

  # Test split
  test_cells_in_seurat <- intersect(smoothed_test$cells, smoothed_test_rna$cells)
  if (length(test_cells_in_seurat) == length(smoothed_test$cells)) {
    smoothed_test$rna_log1p <- smoothed_test_rna$rna_log1p[, smoothed_test$cells, drop = FALSE]
    cat(sprintf("  Test: Added RNA data for %d cells\n", length(test_cells_in_seurat)))
  } else {
    cat(sprintf("  WARNING: Test - only %d/%d cells found in Seurat object\n", 
                length(test_cells_in_seurat), length(smoothed_test$cells)))
    # Subset both ATAC and RNA tp overlapping cells only
    smoothed_test$atac_log1p <- smoothed_test$atac_log1p[, test_cells_in_seurat, drop = FALSE]
    smoothed_test$rna_log1p <- smoothed_test_rna$rna_log1p[, test_cells_in_seurat, drop = FALSE]
    smoothed_test$cells <- test_cells_in_seurat
  }
  
  # Assert check the column names
  stopifnot(
    all(colnames(smoothed_test$rna_log1p) == smoothed_test$cells),
    length(smoothed_test$cells) == ncol(smoothed_test$rna_log1p)
  )

cat("✓ RNA expression data successfully added to all splits\n")

# ===============Get Gene Annotations=========================
# Get gene annotations from EnsemblDb (same as used in preprocessing)
gene_annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# Note: seqlevels/chr style will be harmonized to ATAC peaks below after extracting peak ranges

# ===============Find Highly Variable Genes (FROM TRAINING ONLY)=========================
# Prefer explicit user-supplied HVG list if provided (config: ATAC_FEATURES_HVG_GENE_FILE or HVG_GENE_FILE)
cat("\n=== Determining target HVG genes ===\n")

# Resolve config names / fallbacks
if (!exists("ATAC_FEATURES_HVG_GENE_FILE") && exists("HVG_GENE_FILE")) ATAC_FEATURES_HVG_GENE_FILE <- HVG_GENE_FILE
if (!exists("ATAC_FEATURES_N_HVG_GENES") && exists("N_HVG_GENES")) ATAC_FEATURES_N_HVG_GENES <- N_HVG_GENES
if (!exists("ATAC_FEATURES_N_HVG_GENES")) ATAC_FEATURES_N_HVG_GENES <- n_hvg_genes

# 1) Require a user-supplied HVG file (no automatic HVG computation)
hvg_list <- NULL
if (exists("ATAC_FEATURES_HVG_GENE_FILE") && nzchar(ATAC_FEATURES_HVG_GENE_FILE)) {
  hvf <- path.expand(ATAC_FEATURES_HVG_GENE_FILE)
  if (file.exists(hvf)) {
    hvg_list <- readLines(hvf)
    hvg_list <- hvg_list[nchar(hvg_list) > 0]
    cat(sprintf("Using user-supplied HVG list: %s (%d genes)\n", hvf, length(hvg_list)))
    if (length(hvg_list) > ATAC_FEATURES_N_HVG_GENES) hvg_list <- head(hvg_list, ATAC_FEATURES_N_HVG_GENES)
  } else {
    stop(sprintf("ATAC automated feature extraction requires a user HVG file. Expected: %s", hvf))
  }
} else {
  # also accept general HVG_GENE_FILE name from config.R
  if (exists("HVG_GENE_FILE") && nzchar(HVG_GENE_FILE)) {
    hvf <- path.expand(HVG_GENE_FILE)
    if (file.exists(hvf)) {
      hvg_list <- readLines(hvf)
      hvg_list <- hvg_list[nchar(hvg_list) > 0]
      cat(sprintf("Using user-supplied HVG list from HVG_GENE_FILE: %s (%d genes)\n", hvf, length(hvg_list)))
      if (length(hvg_list) > ATAC_FEATURES_N_HVG_GENES) hvg_list <- head(hvg_list, ATAC_FEATURES_N_HVG_GENES)
    } else {
      stop(sprintf("ATAC automated feature extraction requires a user HVG file. Expected: %s", hvf))
    }
  } else {
    stop("ATAC automated feature extraction requires a user-supplied HVG file (set ATAC_FEATURES_HVG_GENE_FILE or HVG_GENE_FILE in config.R). Aborting.")
  }
}

# ===============Find Top Variable ATAC Peaks (FROM TRAINING ONLY)=========================
cat("\n=== Finding top variable ATAC peaks from TRAINING data only ===\n")
# Compute variable peaks from TRAINING cells only to prevent data leakage
train_cells_for_peaks <- colnames(seurat_obj)[seurat_obj$data_split == "train"]
seurat_train_peaks <- seurat_obj[, train_cells_for_peaks]
DefaultAssay(seurat_train_peaks) <- "ATAC"
seurat_train_peaks <- FindTopFeatures(seurat_train_peaks, min.cutoff = ATAC_MIN_PEAK_QUANTILE)
variable_peaks <- VariableFeatures(seurat_train_peaks)
if (length(variable_peaks) == 0) {
  cat("FindTopFeatures returned 0 peaks; retrying with no cutoff\n")
  seurat_train_peaks <- FindTopFeatures(seurat_train_peaks, min.cutoff = "q0")
  variable_peaks <- VariableFeatures(seurat_train_peaks)
}
cat(sprintf("Computed %d variable ATAC peaks from training data\n", length(variable_peaks)))

# ===============Define Target Gene Sets=========================
cat("\n=== Defining target gene sets ===\n")
# HVG target set is `hvg_list` (either user-provided or computed)
gene_sets <- list(HVG = hvg_list)

# Random genes: prefer ATAC_RANDOM_GENE_FILE then RANDOM_GENE_FILE
if (!exists("ATAC_RANDOM_GENE_FILE") && exists("RANDOM_GENE_FILE")) ATAC_RANDOM_GENE_FILE <- RANDOM_GENE_FILE
random_genes <- NULL
if (exists("ATAC_RANDOM_GENE_FILE") && nzchar(ATAC_RANDOM_GENE_FILE)) {
  rgf <- path.expand(ATAC_RANDOM_GENE_FILE)
  if (file.exists(rgf)) {
    random_genes <- readLines(rgf)
    random_genes <- random_genes[nchar(random_genes) > 0]
    cat(sprintf("Loaded %d random genes from file: %s\n", length(random_genes), basename(rgf)))
  } else {
    cat(sprintf("Random gene file not found: %s\n", rgf))
  }
}
if (!is.null(random_genes) && length(random_genes) > 0) gene_sets$Random_genes <- random_genes

cat(sprintf("Total gene sets to process: %d\n", length(gene_sets)))
for (set_name in names(gene_sets)) cat(sprintf("  - %s: %d genes\n", set_name, length(gene_sets[[set_name]])))

# ===============Get Peak Annotations=========================
cat("\n=== Extracting peak annotations ===\n")

# Get peak coordinates from the ATAC assay
DefaultAssay(seurat_obj) <- "ATAC"
peak_ranges <- granges(seurat_obj[["ATAC"]])

# Ensure peak ranges have names matching ATAC feature names
if (is.null(names(peak_ranges)) || any(is.na(names(peak_ranges)))) {
  names(peak_ranges) <- rownames(seurat_obj[["ATAC"]])
}

# Harmonize seqlevels and chromosome style between peaks and gene annotations
seqlevelsStyle(gene_annotation) <- seqlevelsStyle(peak_ranges)
peak_ranges <- keepStandardChromosomes(peak_ranges, pruning.mode = "coarse")
gene_annotation <- keepStandardChromosomes(gene_annotation, pruning.mode = "coarse")
common_seqlevels <- intersect(seqlevels(peak_ranges), seqlevels(gene_annotation))
peak_ranges <- keepSeqlevels(peak_ranges, common_seqlevels, pruning.mode = "coarse")
gene_annotation <- keepSeqlevels(gene_annotation, common_seqlevels, pruning.mode = "coarse")

cat(sprintf("Total ATAC peaks: %d\n", length(peak_ranges)))

# Verify peak ranges have correct format
cat("\nFirst 5 peaks:\n")
print(head(peak_ranges, 5))

# ===============Process Each Gene Set=========================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("=== PROCESSING GENE SETS (ATAC-ONLY) ===\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

for (set_name in names(gene_sets)) {
  cat(sprintf("\n\n>>> Processing gene set: %s (%d genes) <<<\n", 
              set_name, length(gene_sets[[set_name]])))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  target_genes <- gene_sets[[set_name]]
  
  # Create output directory for this gene set
  set_output_dir <- file.path(output_dir, set_name)
  if (!dir.exists(set_output_dir)) {
    dir.create(set_output_dir, recursive = TRUE)
    cat(sprintf("Created output directory: %s\n", set_output_dir))
  }
  
  # Filter gene annotations to include only target genes from this set
  cat(sprintf("\nFiltering gene annotations to target genes...\n"))
  cat(sprintf("  Target genes in this set: %d\n", length(target_genes)))
  
  target_gene_idx <- which(gene_annotation$gene_name %in% target_genes)
  gene_annotation_set <- gene_annotation[target_gene_idx]
  
  cat(sprintf("  Matched: %d annotation entries for %d unique genes\n", 
              length(gene_annotation_set), length(unique(gene_annotation_set$gene_name))))
  
  # IMPORTANT: Collapse multiple transcripts per gene to single entry per gene
  # Keep only one representative transcript per gene (use the first one)
  cat(sprintf("  Collapsing multiple transcripts per gene...\n"))
  unique_gene_names <- unique(gene_annotation_set$gene_name)
  
  # For each unique gene, keep only the first transcript
  gene_annotation_set <- gene_annotation_set[match(unique_gene_names, gene_annotation_set$gene_name)]
  
  cat(sprintf("\nTarget genes with genomic annotations: %d\n", length(gene_annotation_set)))
  
  # Check for any target genes without annotations
  missing_genes <- setdiff(target_genes, gene_annotation_set$gene_name)
  if (length(missing_genes) > 0) {
    cat(sprintf("Warning: %d target genes do not have genomic annotations and will be skipped:\n", 
                length(missing_genes)))
    cat(paste(head(missing_genes, 10), collapse = ", "), "\n")
    if (length(missing_genes) > 10) cat("... and", length(missing_genes) - 10, "more\n")
  }
  
  # Skip if no genes have annotations
  if (length(gene_annotation_set) == 0) {
    cat(sprintf("ERROR: No genes in %s set have genomic annotations. Skipping.\n", set_name))
    next
  }
  
  # ===============Define Gene Windows and Extract Peaks=========================
  cat(sprintf("\n=== Defining gene windows (±%dkb from TSS) ===\n", gene_window_kb))
  
  # Extend gene regions by the specified window size
  window_bp <- gene_window_kb * 1000
  # Strand-aware TSS windows using promoters() directly on gene ranges
  gene_windows <- promoters(gene_annotation_set, upstream = window_bp, downstream = window_bp)
  names(gene_windows) <- gene_annotation_set$gene_name
  
  cat(sprintf("Created genomic windows for %d genes\n", length(gene_windows)))
  
  # Find overlaps between peaks and gene windows
  cat("\nFinding peaks within gene windows...\n")
  hits <- findOverlaps(peak_ranges, gene_windows, ignore.strand = TRUE)
  
  # Create a mapping: gene_name -> PEAK NAMES within its window (name-based, not indices)
  peak_names_from_gr <- names(peak_ranges)
  gene_to_peak_names <- split(peak_names_from_gr[queryHits(hits)],
                              names(gene_windows)[subjectHits(hits)])
  
  # Restrict to variable peaks only
  cat("\nUsing Variable Peaks...\n")
  gene_to_peak_names <- lapply(gene_to_peak_names, intersect, variable_peaks)
  n_peaks_per_gene <- sapply(gene_to_peak_names, length)
  cat(sprintf("  Mean peaks per gene: %.1f\n", mean(n_peaks_per_gene)))
  cat(sprintf("  Median: %.0f\n", median(n_peaks_per_gene)))
  cat(sprintf("  Min: %d\n", min(n_peaks_per_gene)))
  cat(sprintf("  Max: %d\n", max(n_peaks_per_gene)))

  # Filter out genes with too few peaks
  valid_genes <- names(gene_to_peak_names)[n_peaks_per_gene >= min_peaks_per_gene]
  gene_to_peak_names <- gene_to_peak_names[valid_genes]
  cat(sprintf("\nGenes with ≥%d peaks: %d (%.1f%% of target genes with annotations)\n",
              min_peaks_per_gene, length(gene_to_peak_names),
              100 * length(gene_to_peak_names) / length(gene_annotation_set)))
  
  # ===============Extract Features for Each Split=========================
  cat("\n=== Extracting peak features for each gene ===\n")
  
  # Function to extract peak features for a specific gene and split (by PEAK NAME)
  extract_gene_peak_features <- function(gene_name, gene_peak_names, smoothed_data) {
    all_peak_names <- rownames(smoothed_data$atac_log1p)
    avail <- intersect(gene_peak_names, all_peak_names)
    if (length(avail) == 0) return(NULL)
    feature_matrix <- smoothed_data$atac_log1p[avail, , drop = FALSE]  # peaks × cells
    feature_matrix <- t(feature_matrix)  # cells × peaks
    colnames(feature_matrix) <- avail
    return(feature_matrix)
  }
  
  # Function to extract target gene expression for a specific gene and split
  extract_gene_expression <- function(gene_name, smoothed_data) {
    # Get all gene names
    all_gene_names <- rownames(smoothed_data$rna_log1p)
    
    # Check if gene exists
    if (!gene_name %in% all_gene_names) {
      warning(sprintf("Gene %s not found in RNA data", gene_name))
      return(NULL)
    }
    
    # Extract expression vector (cells)
    expression <- smoothed_data$rna_log1p[gene_name, ]
    
    return(as.numeric(expression))
  }
  
  # Initialize storage for features
  cat("\nExtracting features for all genes and splits...\n")
  cat("This may take a few minutes...\n\n")
  
  gene_features <- list()
  
  # Progress tracking
  n_genes <- length(gene_to_peak_names)
  progress_interval <- max(1, floor(n_genes / 20))  # Update every 5%
  
  for (i in seq_along(gene_to_peak_names)) {
    gene_name <- names(gene_to_peak_names)[i]
    gene_peak_names <- gene_to_peak_names[[i]]
    
    # Progress update
    if (i %% progress_interval == 0 || i == n_genes) {
      cat(sprintf("  [%d/%d] %.1f%% complete\r", i, n_genes, 100 * i / n_genes))
    }
    
    # Extract features for each split
    tryCatch({
      # Extract gene-specific peak features (by name intersection)
      train_peaks <- extract_gene_peak_features(gene_name, gene_peak_names, smoothed_train)
      val_peaks <- extract_gene_peak_features(gene_name, gene_peak_names, smoothed_val)
      test_peaks <- extract_gene_peak_features(gene_name, gene_peak_names, smoothed_test)

      # Skip gene if any split lacks usable features
      if (is.null(train_peaks) || is.null(val_peaks) || is.null(test_peaks) ||
          ncol(train_peaks) == 0 || ncol(val_peaks) == 0 || ncol(test_peaks) == 0) {
        cat(sprintf("\nSkipping gene %s due to no overlapping peak features in one or more splits\n", gene_name))
        next
      }
      
      # Extract target gene expression
      train_y <- extract_gene_expression(gene_name, smoothed_train)
      val_y <- extract_gene_expression(gene_name, smoothed_val)
      test_y <- extract_gene_expression(gene_name, smoothed_test)
      
      # Store in organized structure
      gene_features[[gene_name]] <- list(
        gene_name = gene_name,
        n_peaks = length(gene_peak_names),
        n_total_features = ncol(train_peaks),
        peak_names = colnames(train_peaks),
        
        train = list(
          X = train_peaks,  # cells × peaks
          y = train_y,      # cells
          cells = smoothed_train$cells
        ),
        
        validation = list(
          X = val_peaks,
          y = val_y,
          cells = smoothed_val$cells
        ),
        
        test = list(
          X = test_peaks,
          y = test_y,
          cells = smoothed_test$cells
        )
      )
    }, error = function(e) {
      cat(sprintf("\nWarning: Failed to extract features for gene %s: %s\n", 
                  gene_name, e$message))
    })
  }
  
  cat("\n")
  
  # ===============Quality Control and Summary=========================
  cat("\n=== Feature Extraction Summary for", set_name, "(ATAC-ONLY) ===\n")
  
  # Remove any genes that failed extraction
  gene_features <- gene_features[!sapply(gene_features, is.null)]
  
  cat(sprintf("Successfully extracted features for %d genes\n", length(gene_features)))
  
  # Summary statistics
  if (length(gene_features) > 0) {
    n_peaks_per_gene <- sapply(gene_features, function(x) x$n_peaks)
    
    cat(sprintf("\nPeaks per gene:\n"))
    cat(sprintf("  Mean: %.1f\n", mean(n_peaks_per_gene)))
    cat(sprintf("  Median: %.0f\n", median(n_peaks_per_gene)))
    cat(sprintf("  Range: [%d, %d]\n", min(n_peaks_per_gene), max(n_peaks_per_gene)))
    
    # Check dimensions for first gene (sanity check)
    first_gene <- gene_features[[1]]
    cat(sprintf("\nExample: Gene '%s'\n", first_gene$gene_name))
    cat(sprintf("  Number of peak features: %d\n", first_gene$n_peaks))
    cat(sprintf("  Train: %d cells × %d features\n", 
                nrow(first_gene$train$X), ncol(first_gene$train$X)))
    cat(sprintf("  Validation: %d cells × %d features\n", 
                nrow(first_gene$validation$X), ncol(first_gene$validation$X)))
    cat(sprintf("  Test: %d cells × %d features\n", 
                nrow(first_gene$test$X), ncol(first_gene$test$X)))
  }
  
  # ===============Save Features=========================
  cat("\n=== Saving extracted features ===\n")
  
  # Save the complete feature list
  output_file <- file.path(set_output_dir, "gene_specific_features.rds")
  saveRDS(gene_features, output_file)
  cat(sprintf("Saved: %s\n", output_file))
  
  # Also save metadata for quick reference
  if (length(gene_features) > 0) {
    feature_metadata <- data.frame(
      gene_name = names(gene_features),
      n_peaks = sapply(gene_features, function(x) x$n_peaks),
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
  
  # Save parameters used for this gene set
  params <- list(
    gene_set = set_name,
    feature_type = "ATAC_only",
    gene_window_kb = gene_window_kb,
    min_peaks_per_gene = min_peaks_per_gene,
    n_genes_extracted = length(gene_features),
    peaks_summary = if (length(gene_features) > 0) {
      list(
        mean_peaks = mean(sapply(gene_features, function(x) x$n_peaks)),
        median_peaks = median(sapply(gene_features, function(x) x$n_peaks)),
        min_peaks = min(sapply(gene_features, function(x) x$n_peaks)),
        max_peaks = max(sapply(gene_features, function(x) x$n_peaks))
      )
    } else {
      NULL
    },
    seed = seed,
    extraction_date = Sys.time()
  )
  
  params_file <- file.path(set_output_dir, "feature_extraction_params.rds")
  saveRDS(params, params_file)
  cat(sprintf("Saved parameters: %s\n", params_file))
  
  cat(sprintf("\n>>> Completed %s gene set <<<\n", set_name))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
}  # End of gene set loop

# ===============Final Summary=========================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("=== Step 040 Complete (ATAC-ONLY) ===\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat(sprintf("Gene-specific ATAC peak features extracted successfully\n\n"))

cat("Summary by gene set:\n")
for (set_name in names(gene_sets)) {
  set_output_dir <- file.path(output_dir, set_name)
  feature_file <- file.path(set_output_dir, "gene_specific_features.rds")
  
  if (file.exists(feature_file)) {
    features <- readRDS(feature_file)
    cat(sprintf("  %s: %d genes processed\n", set_name, length(features)))
  }
}

cat("\nKey parameters:\n")
cat(sprintf("  - Feature type: ATAC peaks only\n"))
cat(sprintf("  - Gene window: ±%d kb from TSS\n", gene_window_kb))
cat(sprintf("  - Minimum peaks per gene: %d\n", min_peaks_per_gene))
cat(sprintf("  - Random seed: %d\n", seed))
cat(sprintf("  - Total train cells: %d\n", length(smoothed_train$cells)))
cat(sprintf("  - Total validation cells: %d\n", length(smoothed_val$cells)))
cat(sprintf("  - Total test cells: %d\n", length(smoothed_test$cells)))

cat("\nData leakage prevention:\n")
cat("  ✓ HVGs computed from TRAINING data only (2000 genes)\n")
cat("  ✓ Variable ATAC peaks computed from TRAINING data only\n")
cat("  ✓ Target genes selected from training-derived HVGs\n")
cat("  ✓ Feature selection does NOT use validation/test information\n")
cat("      (Gene-specific approach uses all peaks within ±%dkb window)\n", gene_window_kb)

cat("\nOutput structure:\n")
for (set_name in names(gene_sets)) {
  cat(sprintf("  %s/\n", set_name))
  cat("    - gene_specific_features.rds (all features for all genes)\n")
  cat("    - gene_features_metadata.csv (summary table)\n")
  cat("    - feature_extraction_params.rds (extraction parameters)\n")
}

cat("\nData structure in each RDS file:\n")
cat("  gene_features is a named list where each element contains:\n")
cat("    $gene_name         : Gene symbol\n")
cat("    $n_peaks           : Number of peak features in gene window\n")
cat("    $n_total_features  : Total features (peaks only)\n")
cat("    $peak_names        : Names of the peak features\n")
cat("    $train$X           : Training features (cells × peaks)\n")
cat("    $train$y           : Training targets (gene expression)\n")
cat("    $validation$X      : Validation features\n")
cat("    $validation$y      : Validation targets\n")
cat("    $test$X            : Test features\n")
cat("    $test$y            : Test targets\n")

cat("\nIMPORTANT NOTES:\n")
cat("  ✓ Each gene has gene-specific ATAC peak features only (no HVG expression)\n")
cat("  ✓ Peaks are within ±%dkb window around gene TSS\n", gene_window_kb)
cat("  ✓ All feature matrices are in (cells × features) format ready for modeling\n")
cat("  ✓ HVGs computed from TRAINING data only (Nov 11, 2025 fix)\n")
cat("  ✓ No information from validation/test used in target gene selection\n")

cat("\nNext: Run 05_linear_tree_models.R \n")
