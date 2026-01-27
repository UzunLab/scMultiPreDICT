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

# Target gene directory (from combined pipeline)
INPUT_TARGET_GENES_DIR <- file.path(COMBINED_PIPELINE_DIR, "target_genes", SAMPLE_NAME)

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
# Use the same target genes as the combined pipeline for fair comparison
HVG_GENE_FILE <- file.path(INPUT_TARGET_GENES_DIR, "target_genes_hvg_100.txt")
RANDOM_GENE_FILE <- file.path(INPUT_TARGET_GENES_DIR, "target_genes_random_100.txt")

# Which gene sets to analyze
MODEL_GENE_SET <- "both"  # Options: "HVG", "Random_genes", "both"

# ============================================================================
# SECTION 8: FEATURE EXTRACTION PARAMETERS
# ============================================================================
# For RNA-only, features are top HVG expression values
N_HVG_FEATURES <- 50  # Number of HVG features per target gene

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
# SECTION 10: COMPUTATIONAL RESOURCES  & SLURM SETTINGS
# ============================================================================
N_CORES <- 4
MAX_CORES_TRAINING <- 32
USE_DISK_BACKED <- FALSE

CONDA_ENV_R <- "r-bioc-43"
CONDA_ENV_PYTHON <- "python_env"

SLURM_PARTITION <- "compute"

SLURM_RESOURCES <- list(
  step_030 = list(cpus = 32, mem_gb = 128, time_hours = 6, job_name = "metacell"),
  step_040 = list(cpus = 32, mem_gb = 128, time_hours = 8, job_name = "feature_extract"),
  step_050 = list(cpus = 32, mem_gb = 200, time_hours = 12, job_name = "model_linear"),
  step_060 = list(cpus = 16, mem_gb = 200, time_hours = 24, job_name = "model_nn")
)

SLURM_LOG_DIR <- "~/scRNA_only/jobs/logs"

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
OUTPUT_MODELS_NN_DIR <- file.path(BASE_RESULTS_DIR, "models/NEURAL_NETWORKS", SAMPLE_NAME, DIMRED_METHOD_SUFFIX, NN_ARCH_LABEL)
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

#' Generate a SLURM sbatch file for a pipeline step
#' @param step_name Name of step (e.g., "step_031")
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
#' Generate all sbatch files for the pipeline (scRNA-only)
#' @param output_dir Directory to save sbatch files
generate_all_sbatch <- function(output_dir = ".") {
  cat("\n=== Generating SLURM sbatch files (scRNA-only) ===\n\n")

  # Step 030: Metacell Creation (linear method)
  generate_sbatch("step_030", "03_metacell_creation.R", output_dir)

  # Step 040: Feature Extraction
  generate_sbatch("step_040", "04_feature_extraction.R", output_dir)

  # Step 050: Linear/RF Models
  generate_sbatch("step_050", "05_linear_tree_models.R", output_dir)

  # Step 060: Neural Network
  generate_sbatch("step_060", "06_neural_network.R", output_dir)

  cat("\nAll sbatch files generated in:", output_dir, "\n")
  cat("\n=== Pipeline Submission Order ===\n")
  cat("# Submit in order, waiting for each to complete:\n")
  cat(sprintf("sbatch step_030_%s.sbatch\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_040_%s.sbatch\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_050_%s.sbatch\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_060_%s.sbatch\n", SAMPLE_NAME))
  cat("\n")
}