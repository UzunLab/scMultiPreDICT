# ============================================================================
# RNA-ONLY PIPELINE CONFIGURATION TEMPLATE
# ============================================================================
# This is a TEMPLATE file. To use:
#   1. Copy this file to config.R:  cp config_template.R config.R
#   2. Edit config.R with YOUR dataset-specific paths
#   3. Run the pipeline: Rscript run_pipeline.R
#
# IMPORTANT: The RNA-only pipeline requires outputs from the combined pipeline:
#   - Seurat object with train/val/test splits
#   - Target gene lists (HVG and/or random genes)
#
# Run combined pipeline Steps 1-2 first before using this pipeline.
# ============================================================================

# ============================================================================
# SECTION 1: SAMPLE IDENTIFICATION
# ============================================================================
SAMPLE_NAME <- "YOUR_SAMPLE_NAME"  # Must match the combined pipeline sample name
PROJECT_NAME <- "YOUR_PROJECT_NAME_RNA_only"

# ============================================================================
# SECTION 2: INPUT DATA PATHS (FROM COMBINED PIPELINE)
# ============================================================================
# Path to outputs from the combined pipeline
COMBINED_PIPELINE_DIR <- "/path/to/combined/pipeline/output"

# Seurat object with splits (from combined pipeline Step 02)
INPUT_SEURAT_SPLITS <- file.path(COMBINED_PIPELINE_DIR, "splits", SAMPLE_NAME,
                                  paste0(SAMPLE_NAME, "_seurat_obj_with_splits.rds"))

# ============================================================================
# SECTION 3: OUTPUT DIRECTORIES
# ============================================================================
BASE_OUTPUT_DIR <- "~/scMultiPreDICT_output/rna_only/processed/"
OUTPUT_METACELLS_DIR <- file.path(BASE_OUTPUT_DIR, "metacells", SAMPLE_NAME)
OUTPUT_FEATURES_DIR <- file.path(BASE_OUTPUT_DIR, "features_extracted", SAMPLE_NAME)

# ============================================================================
# SECTION 4: SPECIES/GENOME CONFIGURATION
# ============================================================================
SPECIES <- "mouse"  # or "human"

if (SPECIES == "mouse") {
  GENOME <- "mm10"
  ENSDB_PACKAGE <- "EnsDb.Mmusculus.v79"
}
if (SPECIES == "human") {
  GENOME <- "hg38"
  ENSDB_PACKAGE <- "EnsDb.Hsapiens.v86"
}

# ============================================================================
# SECTION 5: DIMENSIONALITY REDUCTION METHOD
# ============================================================================
DIM_REDUCTION_METHOD <- "PCA"

DIMRED_METHOD_SUFFIX <- switch(
  DIM_REDUCTION_METHOD,
  "PCA" = "pca",
  tolower(DIM_REDUCTION_METHOD)
)

# Update output dirs with method suffix
OUTPUT_METACELLS_DIR <- file.path(BASE_OUTPUT_DIR, "metacells", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
OUTPUT_FEATURES_DIR <- file.path(BASE_OUTPUT_DIR, "features_extracted", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)

# ============================================================================
# SECTION 6: METACELL/KNN SMOOTHING PARAMETERS
# ============================================================================
PCA_DIMS <- 30
N_VARIABLE_FEATURES <- 3000
K_NEIGHBORS <- 20
SEED_METACELL <- 2025

# ============================================================================
# SECTION 7: TARGET GENE CONFIGURATION
# ============================================================================
# Use the SAME pre-computed target gene lists across all 3 pipelines for fair comparison.
# Pre-computed files are provided in the repo under data/target_genes/.
# Use ABSOLUTE paths (or paths relative to the directory you run Rscript from).
HVG_GENE_FILE <- "<PATH_TO_HVG_GENE_LIST>"       # e.g., "/abs/path/to/repo/data/target_genes/E7.5_rep2/target_genes_hvg_100.txt"
RANDOM_GENE_FILE <- "<PATH_TO_RANDOM_GENE_LIST>"  # e.g., "/abs/path/to/repo/data/target_genes/E7.5_rep2/target_genes_random_100.txt"

# Which gene sets to analyze
MODEL_GENE_SET <- "both"  # Options: "HVG", "Random_genes", "both"

# ============================================================================
# SECTION 8: FEATURE EXTRACTION PARAMETERS
# ============================================================================
# For RNA-only, features are top HVG expression values
N_HVG_FEATURES <- 50  # Number of HVG features per target gene
SEED_FEATURES <- 123  # Random seed for feature extraction

# ============================================================================
# SECTION 9: MODEL TRAINING PARAMETERS
# ============================================================================
SEED_MODEL <- 123
RF_N_TREES <- 500

# Neural network settings
NN_HIDDEN_UNITS <- 256
NN_N_HIDDEN_LAYERS <- 3
NN_DROPOUT_RATE <- 0.1
NN_LEARNING_RATE <- 0.001
NN_BATCH_SIZE <- 256
NN_MAX_EPOCHS <- 100
NN_EARLY_STOP_PATIENCE <- 10

# Grid search options
NN_USE_GRID_SEARCH <- TRUE
NN_GRID_UNITS <- c(256)
NN_GRID_DROPOUT <- c(0, 0.1)
NN_GRID_BATCH <- c(256, 512)

# Conda environment for TensorFlow/Keras
CONDA_ENV_NAME <- "your_conda_env"
# ============================================================================
# SECTION 10: COMPUTATIONAL RESOURCES
# ============================================================================
N_CORES <- 4
MAX_CORES_TRAINING <- 32
USE_DISK_BACKED <- FALSE

# ============================================================================
# SECTION 11: RESULTS OUTPUT DIRECTORIES
# ============================================================================
BASE_RESULTS_DIR <- "~/scMultiPreDICT_output/rna_only/results/"

NN_ARCH_LABEL <- switch(
  as.character(NN_N_HIDDEN_LAYERS),
  "1" = "One_hidden_layer",
  "2" = "Two_hidden_layer", 
  "3" = "Three_hidden_layer",
  paste0(NN_N_HIDDEN_LAYERS, "_hidden_layer")
)

OUTPUT_MODELS_LINEAR_DIR <- file.path(BASE_RESULTS_DIR, "models/LINEAR_AND_TREE_BASED", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
OUTPUT_MODELS_NN_DIR <- file.path(BASE_RESULTS_DIR, "models/NEURAL_NETWORKS", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
OUTPUT_FIGURES_DIR <- file.path(BASE_RESULTS_DIR, "figures", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

create_output_directories <- function() {
  dirs <- c(
    OUTPUT_METACELLS_DIR,
    OUTPUT_FEATURES_DIR,
    OUTPUT_MODELS_LINEAR_DIR,
    OUTPUT_MODELS_NN_DIR,
    OUTPUT_FIGURES_DIR
  )
  for (dir in dirs) {
    dir <- path.expand(dir)
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat(sprintf("Created directory: %s\n", dir))
    }
  }
}

load_annotation_packages <- function() {
  if (SPECIES == "mouse") {
    suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
  } else if (SPECIES == "human") {
    suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
  }
}

print_config <- function() {
  cat("\n=======================================================================\n")
  cat("RNA-ONLY PIPELINE CONFIGURATION\n")
  cat("=======================================================================\n\n")
  cat("Sample:", SAMPLE_NAME, "\n")
  cat("Species:", SPECIES, "(", GENOME, ")\n")
  cat("Dimensionality Reduction:", DIM_REDUCTION_METHOD, "\n")
  cat("Input Seurat:", INPUT_SEURAT_SPLITS, "\n\n")
}

print_output_directories <- function() {
  cat("\n=======================================================================\n")
  cat("OUTPUT DIRECTORIES\n")
  cat("=======================================================================\n\n")
  cat("  Metacells:       ", path.expand(OUTPUT_METACELLS_DIR), "\n")
  cat("  Features:        ", path.expand(OUTPUT_FEATURES_DIR), "\n")
  cat("  Linear Models:   ", path.expand(OUTPUT_MODELS_LINEAR_DIR), "\n")
  cat("  Neural Networks: ", path.expand(OUTPUT_MODELS_NN_DIR), "\n")
  cat("  Figures:         ", path.expand(OUTPUT_FIGURES_DIR), "\n\n")
}

validate_config <- function() {
  errors <- c()
  if (!file.exists(path.expand(INPUT_SEURAT_SPLITS))) {
    errors <- c(errors, paste("INPUT_SEURAT_SPLITS not found:", INPUT_SEURAT_SPLITS))
  }
  if (length(errors) > 0) {
    cat("Configuration errors:\n")
    for (err in errors) cat("  ERROR:", err, "\n")
    stop("Please fix configuration errors before proceeding.")
  }
  cat("Configuration validated successfully!\n")
}

cat("\n[CONFIG] Loaded RNA-only configuration for sample:", SAMPLE_NAME, "\n")
cat("[CONFIG] Species:", SPECIES, "| Genome:", GENOME, "\n\n")


# Flag to prevent re-sourcing when called via run_pipeline.R
CONFIG_LOADED <- TRUE
