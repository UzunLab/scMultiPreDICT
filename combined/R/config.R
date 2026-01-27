# ============================================================================
# PIPELINE CONFIGURATION TEMPLATE
# ============================================================================
# This is a TEMPLATE file. To use:
#   1. Copy this file to config.R:  cp config_template.R config.R
#   2. Edit config.R with YOUR dataset-specific paths and parameters
#   3. Run the pipeline: Rscript run_pipeline.R
#
# DO NOT edit this template directly - keep it as a reference.
# ============================================================================

# ============================================================================
# SECTION 1: SAMPLE IDENTIFICATION
# ============================================================================

# Sample name - used for file naming and output organization
SAMPLE_NAME <- "<YOUR_SAMPLE_NAME>"  # e.g., "E7.5_rep1", "T_Cells", "PBMC_10k"

# Project name (optional, for documentation)
PROJECT_NAME <- "<YOUR_PROJECT_NAME>"

# ============================================================================
# SECTION 2: INPUT DATA PATHS
# ============================================================================
# Paths to the raw data files (10X Genomics format)
# Use absolute paths or paths relative to your home directory (~/)


# RNA data paths (Matrix Market format)
INPUT_MTX <- "<PATH_TO_MATRIX_MTX>"  # e.g., "/path/to/your/matrix.mtx.gz"
INPUT_FEATURES <- "<PATH_TO_FEATURES_TSV>"  # e.g., "/path/to/your/features.tsv.gz"
INPUT_BARCODES <- "<PATH_TO_BARCODES_TSV>"  # e.g., "/path/to/your/barcodes.tsv.gz"

# ATAC fragments file (must be sorted and indexed)
INPUT_FRAGMENTS <- "<PATH_TO_FRAGMENTS_TSV>"  # e.g., "/path/to/your/fragments.tsv.gz"

# ============================================================================
# SECTION 3: OUTPUT DIRECTORIES
# ============================================================================

# Base output directory - all outputs will be organized under this
BASE_OUTPUT_DIR <- "<PATH_TO_OUTPUT_DIRECTORY>"  # e.g., "~/scMultiPreDICT_output/processed/"

# ============================================================================
# SECTION 3a: DIMENSIONALITY REDUCTION METHOD
# ============================================================================
# Choose which dimensionality reduction approach to use for metacell construction
# Options:
#   "pca_lsi"      - PCA for RNA, LSI for ATAC (default, fast, done in R)
#   "wnn"          - Weighted Nearest Neighbors from Seurat (done in R)
#   "multivi"      - Joint MultiVI autoencoder (requires Python)
#   "scvi_peakvi"  - Separate scVI + PeakVI autoencoders (requires Python)

DIM_REDUCTION_METHOD <- "pca_lsi"

# Short suffix for folder names (auto-derived from DIM_REDUCTION_METHOD)
DIMRED_METHOD_SUFFIX <- switch(
  DIM_REDUCTION_METHOD,
  "pca_lsi" = "pca_lsi",
  "wnn" = "wnn",
  "multivi" = "multivi",
  "scvi_peakvi" = "scvi_peakvi",
  DIM_REDUCTION_METHOD
)

# ----------------------------------------------------------------------------
# Derived output directories (automatically constructed)
# ----------------------------------------------------------------------------
# Steps 01-02: Common preprocessing (NOT method-specific)
OUTPUT_SEURAT_DIR <- file.path(BASE_OUTPUT_DIR, "seurat_obj", SAMPLE_NAME)
OUTPUT_SPLITS_DIR <- file.path(BASE_OUTPUT_DIR, "splits", SAMPLE_NAME)
OUTPUT_TARGET_GENES_DIR <- file.path(BASE_OUTPUT_DIR, "target_genes", SAMPLE_NAME)
OUTPUT_LATENT_DIR <- file.path(BASE_OUTPUT_DIR, "latent_spaces", SAMPLE_NAME)

# Steps 03+: Method-specific outputs
OUTPUT_METACELLS_DIR <- file.path(BASE_OUTPUT_DIR, "metacells", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
OUTPUT_FEATURES_DIR <- file.path(BASE_OUTPUT_DIR, "features_extracted", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)

# ============================================================================
# SECTION 4: SPECIES/GENOME CONFIGURATION
# ============================================================================

# Supported options: "mouse" or "human"
SPECIES <- "<SPECIES>"  # "mouse" or "human"

# Species-specific settings (auto-configured)
if (SPECIES == "mouse") {
  GENOME <- "mm10"
  ENSDB_PACKAGE <- "EnsDb.Mmusculus.v79"
  BSGENOME_PACKAGE <- "BSgenome.Mmusculus.UCSC.mm10"
  MT_PATTERN <- "^mt-"
}

if (SPECIES == "human") {
  GENOME <- "hg38"
  ENSDB_PACKAGE <- "EnsDb.Hsapiens.v86"
  BSGENOME_PACKAGE <- "BSgenome.Hsapiens.UCSC.hg38"
  MT_PATTERN <- "^MT-"
}

# ============================================================================
# SECTION 5: QUALITY CONTROL PARAMETERS
# ============================================================================
# Adjust these thresholds based on your data quality

# RNA QC thresholds
QC_MIN_FEATURES_RNA <- 1000      # Minimum genes detected per cell
QC_MAX_FEATURES_RNA <- 6000      # Maximum genes detected per cell (doublet filter)
QC_MAX_PERCENT_MT <- 35          # Maximum mitochondrial percentage

# ATAC QC thresholds
QC_MIN_COUNT_ATAC <- 1000        # Minimum ATAC fragments per cell
QC_MAX_COUNT_ATAC <- 100000      # Maximum ATAC fragments per cell
QC_MIN_COUNT_RNA <- 1000         # Minimum RNA counts per cell
QC_MAX_COUNT_RNA <- 40000        # Maximum RNA counts per cell
QC_MAX_NUCLEOSOME_SIGNAL <- 1.7  # Maximum nucleosome signal
QC_MIN_TSS_ENRICHMENT <- 3       # Minimum TSS enrichment score

# Gene expression filter: Keep genes expressed in at least this fraction of cells
QC_MIN_GENE_FRACTION <- 0.10

# ============================================================================
# SECTION 6: DATA SPLITTING PARAMETERS
# ============================================================================
# Train/Validation/Test split proportions (must sum to 1.0)
SPLIT_TRAIN <- 0.70
SPLIT_VAL <- 0.20
SPLIT_TEST <- 0.10

# Random seed for reproducibility
SEED_SPLIT <- 42

# ============================================================================
# SECTION 7: METACELL/KNN SMOOTHING PARAMETERS
# ============================================================================
# Dimensionality reduction parameters
PCA_DIMS <- 30                   # Number of PCA dimensions for RNA
LSI_DIMS <- 30                   # Number of LSI dimensions for ATAC
N_VARIABLE_FEATURES <- 3000      # Number of highly variable genes for PCA

# k-NN smoothing parameters
K_NEIGHBORS <- 20                # Number of neighbors for k-NN smoothing

# Random seed for reproducibility
SEED_METACELL <- 2025

# ============================================================================
# SECTION 8: TARGET GENE CONFIGURATION
# ============================================================================
# Target genes define WHICH genes' expression to predict (the response variable Y).
#
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  OPTION 1: Use pre-computed files (REQUIRED for reproducing results)    ║
# ║  OPTION 2: Auto-select (ONLY for analyzing NEW datasets)                ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# ---- OPTION 1: Pre-computed Target Gene Files ----
# For reproducing published results, use the provided files:
#   E7.5_rep1: "data/target_genes/E7.5_rep1/target_genes_hvg_100.txt"
#   E7.5_rep2: "data/target_genes/E7.5_rep2/target_genes_hvg_100.txt"
#   T_Cells:   "data/target_genes/T_Cells/target_genes_hvg_100.txt"
#
# Set to empty string "" to use auto-selection (Option 2)
HVG_GENE_FILE <- ""
RANDOM_GENE_FILE <- ""

# ---- OPTION 2: Auto-selection Parameters (for NEW datasets) ----
N_HVG_GENES <- 100               # Number of HVGs to use as target genes
N_RANDOM_TARGET_GENES <- 100     # Number of random non-HVG genes
TARGET_MIN_DETECTION <- 5        # Min % of cells expressing the gene
TARGET_MAX_DETECTION <- 95       # Max % of cells expressing the gene
TARGET_MIN_MEAN_EXPR <- 0.1      # Min mean expression (log-normalized)
SEED_TARGET_GENES <- 2025        # Random seed for target gene selection

# Path for auto-generated target gene list (created by Step 02b)
TARGET_GENE_FILE <- file.path(OUTPUT_TARGET_GENES_DIR, 
                               sprintf("target_genes_random_%d.txt", N_RANDOM_TARGET_GENES))

SEED_FEATURES <- 123        # Random seed for target gene selection
# ============================================================================
# SECTION 9: AUTOENCODER LATENT SPACE PATHS (for multivi, scvi_peakvi methods)
# ============================================================================
# These paths are auto-constructed based on OUTPUT_LATENT_DIR
MULTIVI_LATENT_PATH <- file.path(OUTPUT_LATENT_DIR, "latent_multivi_all.csv")
SCVI_LATENT_PATH <- file.path(OUTPUT_LATENT_DIR, "latent_scvi_rna_all.csv")
PEAKVI_LATENT_PATH <- file.path(OUTPUT_LATENT_DIR, "latent_peakvi_atac_all.csv")

# ============================================================================
# SECTION 10: FEATURE EXTRACTION PARAMETERS
# ============================================================================
# Genomic window around TSS for ATAC peak extraction
GENE_WINDOW_KB <- 250            # ±250kb window from TSS

# Minimum peaks required per gene
MIN_PEAKS_PER_GENE <- 5

# Number of HVG expression features to include
N_HVG_FEATURES <- 1000  # Number of top HVGs to use as expression features

# ============================================================================
# SECTION 11: MODEL TRAINING PARAMETERS
# ============================================================================
# Which gene sets to train models on
MODEL_GENE_SET <- "both"         # Options: "HVG", "Random_genes", "both"

# Random seed for model training
SEED_MODEL <- 123

# --- Random Forest ---
RF_N_TREES <- 500

# --- Deep Learning (Step 06) ---
NN_HIDDEN_UNITS <- 256           # Units in first hidden layer
NN_N_HIDDEN_LAYERS <- 3          # Number of hidden layers (1, 2, or 3)
NN_DROPOUT_RATE <- 0.1           # Dropout rate for regularization
NN_LEARNING_RATE <- 0.001        # Adam optimizer learning rate
NN_BATCH_SIZE <- 256             # Training batch size
NN_MAX_EPOCHS <- 100             # Maximum training epochs
NN_EARLY_STOP_PATIENCE <- 10     # Early stopping patience

# Grid search options (set to FALSE to use fixed params above)
NN_USE_GRID_SEARCH <- TRUE
NN_GRID_UNITS <- c(256)
NN_GRID_DROPOUT <- c(0, 0.1)
NN_GRID_BATCH <- c(256, 512)


# Conda environment for TensorFlow/Keras
CONDA_ENV_NAME <- "<YOUR_CONDA_ENV>"

# ============================================================================
# SECTION 12: COMPUTATIONAL RESOURCES
# ============================================================================
N_CORES <- 4                     # Default cores for local runs
MAX_CORES_TRAINING <- 32         # Maximum cores for model training
USE_DISK_BACKED <- FALSE         # Use disk-backed matrices for large datasets


# Conda environment names
CONDA_ENV_R <- "<YOUR_R_CONDA_ENV>"
CONDA_ENV_PYTHON <- "<YOUR_PYTHON_CONDA_ENV>"


# SLURM settings
SLURM_PARTITION <- "<YOUR_SLURM_PARTITION>"  # e.g., "compute"
SLURM_LOG_DIR <- "<PATH_TO_SLURM_LOG_DIR>"  # e.g., "~/logs"


# ============================================================================
# SECTION 13: RESULTS OUTPUT DIRECTORIES
# ============================================================================
BASE_RESULTS_DIR <- "<PATH_TO_RESULTS_DIRECTORY>"  # e.g., "~/scMultiPreDICT_output/results/"

# Model output directories
OUTPUT_MODELS_LINEAR_DIR <- file.path(BASE_RESULTS_DIR, "models/LINEAR_AND_TREE_BASED", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)

NN_ARCH_LABEL <- switch(
  as.character(NN_N_HIDDEN_LAYERS),
  "1" = "One_hidden_layer",
  "2" = "Two_hidden_layer",
  "3" = "Three_hidden_layer",
  paste0(NN_N_HIDDEN_LAYERS, "_hidden_layer")
)

OUTPUT_MODELS_NN_DIR <- file.path(BASE_RESULTS_DIR, "models/NEURAL_NETWORKS", SAMPLE_NAME, DIMRED_METHOD_SUFFIX, NN_ARCH_LABEL)
OUTPUT_SHAP_DIR <- file.path(BASE_RESULTS_DIR, "shap_analysis", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
OUTPUT_FIGURES_DIR <- file.path(BASE_RESULTS_DIR, "figures", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)

# ============================================================================
# SLURM RESOURCE REQUIREMENTS (adjust based on your dataset size)
# ============================================================================
SLURM_RESOURCES <- list(
  step_01 = list(cpus = 8, mem_gb = 64, time_hours = 4, job_name = "qc_preprocess"),
  step_02b = list(cpus = 4, mem_gb = 32, time_hours = 1, job_name = "target_genes"),
  step_020 = list(cpus = 4, mem_gb = 32, time_hours = 1, job_name = "data_split"),
  step_025 = list(cpus = 8, mem_gb = 64, time_hours = 6, job_name = "train_autoencoder", gpu = TRUE, gpu_count = 1),
  step_03b = list(cpus = 4, mem_gb = 64, time_hours = 1, job_name = "export_mudata"),
  step_030 = list(cpus = 32, mem_gb = 256, time_hours = 8, job_name = "metacell"),
  step_040 = list(cpus = 32, mem_gb = 256, time_hours = 12, job_name = "feature_extract"),
  step_050 = list(cpus = 32, mem_gb = 128, time_hours = 12, job_name = "model_linear"),
  step_060 = list(cpus = 16, mem_gb = 64, time_hours = 24, job_name = "model_nn")
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Create all output directories
create_output_directories <- function() {
  dirs <- c(
    OUTPUT_SEURAT_DIR,
    OUTPUT_SPLITS_DIR,
    OUTPUT_METACELLS_DIR,
    OUTPUT_FEATURES_DIR,
    OUTPUT_TARGET_GENES_DIR,
    OUTPUT_LATENT_DIR,
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

#' Print all configured output directories
#' Useful for verifying paths before running pipeline
print_output_directories <- function() {
  cat("\n")
  cat("=" , rep("=", 70), "\n", sep = "")
  cat("OUTPUT DIRECTORY CONFIGURATION\n")
  cat("=" , rep("=", 70), "\n\n", sep = "")
  
  cat(sprintf("Sample: %s\n", SAMPLE_NAME))
  cat(sprintf("Dimensionality Reduction Method: %s (%s)\n\n", 
              DIM_REDUCTION_METHOD, DIMRED_METHOD_SUFFIX))
  
  cat("COMMON DIRECTORIES (method-independent):\n")
  cat(sprintf("  Seurat Objects:     %s\n", path.expand(OUTPUT_SEURAT_DIR)))
  cat(sprintf("  Data Splits:        %s\n", path.expand(OUTPUT_SPLITS_DIR)))
  cat(sprintf("  Target Genes:       %s\n", path.expand(OUTPUT_TARGET_GENES_DIR)))
  cat(sprintf("  Latent Spaces:      %s\n\n", path.expand(OUTPUT_LATENT_DIR)))
  
  cat("METHOD-SPECIFIC DIRECTORIES (organized by dimred method):\n")
  cat(sprintf("  Metacells:          %s\n", path.expand(OUTPUT_METACELLS_DIR)))
  cat(sprintf("  Features:           %s\n", path.expand(OUTPUT_FEATURES_DIR)))
  cat(sprintf("  Linear Models:      %s\n", path.expand(OUTPUT_MODELS_LINEAR_DIR)))
  cat(sprintf("  Neural Networks:    %s\n", path.expand(OUTPUT_MODELS_NN_DIR)))
  cat(sprintf("  SHAP Analysis:      %s\n", path.expand(OUTPUT_SHAP_DIR)))
  cat(sprintf("  Figures:            %s\n\n", path.expand(OUTPUT_FIGURES_DIR)))
  
  cat("=" , rep("=", 70), "\n\n", sep = "")
}

#' Load appropriate annotation packages
load_annotation_packages <- function() {
  if (SPECIES == "mouse") {
    suppressPackageStartupMessages({
      library(EnsDb.Mmusculus.v79)
      library(BSgenome.Mmusculus.UCSC.mm10)
    })
  } else if (SPECIES == "human") {
    suppressPackageStartupMessages({
      library(EnsDb.Hsapiens.v86)
      library(BSgenome.Hsapiens.UCSC.hg38)
    })
  }
}

#' Get EnsDb object based on species
get_ensdb <- function() {
  if (SPECIES == "mouse") {
    return(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
  } else if (SPECIES == "human") {
    return(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
  }
}

#' Print configuration summary
print_config <- function() {
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("PIPELINE CONFIGURATION\n")
  cat("=", rep("=", 70), "\n\n", sep = "")
  
  cat("Sample:", SAMPLE_NAME, "\n")
  cat("Species:", SPECIES, "(", GENOME, ")\n")
  cat("Dimensionality Reduction:", DIM_REDUCTION_METHOD, "\n")
  cat("Base output:", BASE_OUTPUT_DIR, "\n\n")
  
  cat("QC Thresholds:\n")
  cat("  RNA features:", QC_MIN_FEATURES_RNA, "-", QC_MAX_FEATURES_RNA, "\n")
  cat("  MT percentage: <", QC_MAX_PERCENT_MT, "%\n")
  cat("  ATAC counts:", QC_MIN_COUNT_ATAC, "-", QC_MAX_COUNT_ATAC, "\n\n")
  
  cat("Data Splits:\n")
  cat("  Train:", SPLIT_TRAIN * 100, "%\n")
  cat("  Validation:", SPLIT_VAL * 100, "%\n")
  cat("  Test:", SPLIT_TEST * 100, "%\n\n")
  
  cat("Metacell Parameters:\n")
  cat("  PCA/LSI dims:", PCA_DIMS, "\n")
  cat("  k-NN neighbors:", K_NEIGHBORS, "\n\n")
}

#' Validate configuration
validate_config <- function() {
  errors <- c()
  
  # Check split proportions
  if (abs(SPLIT_TRAIN + SPLIT_VAL + SPLIT_TEST - 1.0) > 0.001) {
    errors <- c(errors, "Split proportions must sum to 1.0")
  }
  
  # Check species
  if (!SPECIES %in% c("mouse", "human")) {
    errors <- c(errors, "SPECIES must be 'mouse' or 'human'")
  }
  
  # Check input files exist
  if (!file.exists(path.expand(INPUT_MTX))) {
    errors <- c(errors, paste("INPUT_MTX not found:", INPUT_MTX))
  }
  if (!file.exists(path.expand(INPUT_FEATURES))) {
    errors <- c(errors, paste("INPUT_FEATURES not found:", INPUT_FEATURES))
  }
  if (!file.exists(path.expand(INPUT_BARCODES))) {
    errors <- c(errors, paste("INPUT_BARCODES not found:", INPUT_BARCODES))
  }
  if (!file.exists(path.expand(INPUT_FRAGMENTS))) {
    errors <- c(errors, paste("INPUT_FRAGMENTS not found:", INPUT_FRAGMENTS))
  }
  
  if (length(errors) > 0) {
    cat("Configuration errors:\n")
    for (err in errors) {
      cat("  ERROR:", err, "\n")
    }
    stop("Please fix configuration errors before proceeding.")
  }
  
  cat("Configuration validated successfully!\n")
}

#' Generate SLURM sbatch file for a pipeline step
#' @param step_name Name of step (e.g., "010", "020")
#' @param script_name Name of R/Python script to run
#' @param output_dir Directory to save sbatch file (default: current directory)
#' @param conda_env Conda environment name (default: from config)
#' @param extra_args Additional arguments to pass to the script
#' @param array_range Optional SLURM array range (e.g., "1-100")
#' @param is_python Whether the script is Python (default: FALSE for R)
generate_sbatch <- function(step_name, script_name, output_dir = ".", 
                            conda_env = NULL, extra_args = "", 
                            array_range = NULL, is_python = FALSE) {
  
  # Get resource settings
  if (!step_name %in% names(SLURM_RESOURCES)) {
    stop(sprintf("Unknown step: %s. Available: %s", 
                 step_name, paste(names(SLURM_RESOURCES), collapse = ", ")))
  }
  res <- SLURM_RESOURCES[[step_name]]
  
  # Default conda environment
  if (is.null(conda_env)) {
    conda_env <- if (is_python) CONDA_ENV_PYTHON else CONDA_ENV_R
  }
  
  # Format time
  hours <- floor(res$time_hours)
  minutes <- round((res$time_hours - hours) * 60)
  time_str <- sprintf("%02d:%02d:00", hours, minutes)
  
  # Log directory
  log_dir <- file.path(path.expand(SLURM_LOG_DIR), res$job_name)
  
  # Check if GPU is requested
  use_gpu <- !is.null(res$gpu) && res$gpu
  gpu_count <- if (!is.null(res$gpu_count)) res$gpu_count else 1
  
  # Build sbatch content
  sbatch_lines <- c(
    "#!/bin/bash",
    sprintf("#SBATCH --job-name=%s_%s", res$job_name, SAMPLE_NAME),
    sprintf("#SBATCH -p %s", if (use_gpu) "gpu" else SLURM_PARTITION),
    sprintf("#SBATCH --time=%s", time_str),
    sprintf("#SBATCH --output=%s/%%x_%%j.log", log_dir),
    sprintf("#SBATCH --error=%s/%%x_%%j.err", log_dir),
    "#SBATCH -N 1",
    "#SBATCH --ntasks-per-node=1",
    sprintf("#SBATCH -c %d", res$cpus),
    sprintf("#SBATCH --mem=%dG", res$mem_gb)
  )
  
  # Add GPU if requested
  if (use_gpu) {
    sbatch_lines <- c(sbatch_lines, sprintf("#SBATCH --gres=gpu:%d", gpu_count))
  }
  
  # Add array if specified
  if (!is.null(array_range)) {
    sbatch_lines <- c(sbatch_lines, sprintf("#SBATCH --array=%s", array_range))
  }
  
  sbatch_lines <- c(sbatch_lines, "", "set -euo pipefail", "")
  
  # Create log directory
  sbatch_lines <- c(sbatch_lines,
    "# Create log directory",
    sprintf("mkdir -p \"%s\"", log_dir),
    ""
  )
  
  # Conda activation
  sbatch_lines <- c(sbatch_lines,
    "# Initialize and activate conda environment",
    "source ~/miniconda3/etc/profile.d/conda.sh",
    sprintf("conda activate %s", conda_env),
    ""
  )
  
  # Environment info
  sbatch_lines <- c(sbatch_lines,
    "# Print environment info",
    "echo \"============================================\"",
    "echo \"Job Information\"",
    "echo \"============================================\"",
    "echo \"Job ID: $SLURM_JOB_ID\"",
    "echo \"Job Name: $SLURM_JOB_NAME\"",
    sprintf("echo \"CPUs: %d\"", res$cpus),
    sprintf("echo \"Memory: %dG\"", res$mem_gb),
    "echo \"Conda Env: $CONDA_DEFAULT_ENV\"",
    "echo \"Start Time: $(date)\"",
    "echo \"============================================\"",
    ""
  )
  
  # Set R library paths if R script
  if (!is_python) {
    sbatch_lines <- c(sbatch_lines,
      "# Set R library paths",
      "export R_LIBS=\"$CONDA_PREFIX/lib/R/library\"",
      "export R_LIBS_USER=\"$CONDA_PREFIX/lib/R/library\"",
      "export R_LIBS_SITE=\"\"",
      ""
    )
  }
  
  # Run script
  if (is_python) {
    run_cmd <- sprintf("python %s %s", script_name, extra_args)
  } else {
    run_cmd <- sprintf("Rscript %s %s", script_name, extra_args)
  }
  
  sbatch_lines <- c(sbatch_lines,
    "# Run the script",
    run_cmd,
    "",
    "# Print completion",
    "echo \"\"",
    "echo \"Job completed at $(date)\"",
    "echo \"Exit code: $?\""
  )
  
  # Write to file
  sbatch_filename <- file.path(output_dir, sprintf("%s_%s.sbatch", step_name, SAMPLE_NAME))
  writeLines(sbatch_lines, sbatch_filename)
  
  cat(sprintf("Generated: %s\n", sbatch_filename))
  cat(sprintf("  Resources: %d CPUs, %dG RAM, %s time\n", res$cpus, res$mem_gb, time_str))
  
  invisible(sbatch_filename)
}

#' Generate all sbatch files for the pipeline
#' @param output_dir Directory to save sbatch files
#' @param gene_set Which gene set to use: "HVG" or "Random_genes"
#' @param output_dir Directory to write sbatch files to
generate_all_sbatch <- function(output_dir = ".") {
  cat("\n=== Generating SLURM sbatch files ===\n\n")
  
  # Step 01: QC and Preprocessing
  generate_sbatch("step_01", "01_quality_control.R", output_dir)
  
  # Step 02b: Select Target Genes
  generate_sbatch("step_02b", "02b_select_target_genes.R", output_dir)
  
  # Step 020: Data Splitting
  generate_sbatch("step_020", "02_data_splitting.R", output_dir)
  
  # Step 03b: Export Seurat to MuData (required for Python steps)
  generate_sbatch("step_03b", "03b_export_seurat_to_mudata.R", output_dir)
  
  # Step 025: Dimensionality Reduction with Autoencoders (Python - GPU)
  # Build arguments for the Python script
  train_autoencoder_input <- file.path(path.expand(OUTPUT_SPLITS_DIR), paste0(SAMPLE_NAME, "_seurat_obj_with_splits_mu.h5mu"))
  train_autoencoder_output <- file.path(path.expand(OUTPUT_METACELLS_DIR), "Autoencoders_For_Dimensionality_Reduction_Projected")
  train_autoencoder_args <- sprintf("--input '%s' --output-dir '%s' --sample-name '%s' --n-latent 30 --max-epochs 500",
                          train_autoencoder_input, train_autoencoder_output, SAMPLE_NAME)
  generate_sbatch("step_025", "train_autoencoders.py", 
                  output_dir, conda_env = CONDA_ENV_PYTHON, is_python = TRUE, extra_args = train_autoencoder_args)
  
  # Step 030: Metacell Creation (linear method)
  generate_sbatch("step_030", "03_metacell_creation_pca_lsi.R", output_dir)
  
  # Step 030 variants (uncomment the one you want to use):
  # --- Autoencoder latent spaces ---
  # generate_sbatch("step_030_multivi", "03_metacell_creation_multivi.R", output_dir)
  # generate_sbatch("step_030_scvi_peakvi", "03_metacell_creation_scvi_peakvi.R", output_dir)
  # --- WMN--
  # generate_sbatch("step_030_split_independent", "03_metacell_creation_wnn.R", output_dir)
  
  # Step 040: Feature Extraction
  generate_sbatch("step_040", "04_feature_extraction.R", output_dir)
  
  # Step 050: Linear/RF Models (uses parallel processing internally)
  generate_sbatch("step_050", "05_linear_tree_models.R", output_dir)
  
  # Step 060: Neural Network
  generate_sbatch("step_060", "06_neural_network.R", output_dir)
  
  cat("\nAll sbatch files generated in:", output_dir, "\n")
  cat("\n=== Pipeline Submission Order ===\n")
  cat("# Submit in order, waiting for each to complete:\n")
  cat(sprintf("sbatch step_01_%s.sbatch      # QC and Preprocessing\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_02b_%s.sbatch     # Select Target Genes\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_020_%s.sbatch     # Data Splitting\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_03b_%s.sbatch     # Export to MuData for Python\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_025_%s.sbatch     # Optional: GPU autoencoder training\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_030_%s.sbatch     # Metacell Creation\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_040_%s.sbatch     # Feature Extraction\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_050_%s.sbatch     # Linear/Tree Models\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_060_%s.sbatch     # Neural Network\n", SAMPLE_NAME))
  cat("\n")
}

# ============================================================================
# AUTO-LOAD MESSAGE
# ============================================================================
cat("\n[CONFIG] Loaded configuration for sample:", SAMPLE_NAME, "\n")
cat("[CONFIG] Species:", SPECIES, "| Genome:", GENOME, "\n")
cat("[CONFIG] Dimensionality Reduction:", DIM_REDUCTION_METHOD, "\n")
cat("[CONFIG] Use print_config() for full configuration\n\n")
