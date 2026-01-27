# ============================================================================
# ATAC-ONLY PIPELINE CONFIGURATION TEMPLATE
# ============================================================================
# This is a TEMPLATE file. To use:
#   1. Copy this file to config.R:  cp config_template.R config.R
#   2. Edit config.R with YOUR dataset-specific paths
#   3. Run the pipeline: Rscript run_pipeline.R
#
# IMPORTANT: The ATAC-only pipeline requires outputs from the combined pipeline:
#   - Seurat object with train/val/test splits
#   - Target gene lists (HVG and/or random genes)
#
# Run combined pipeline Steps 1-2 first before using this pipeline.
# ============================================================================

# ============================================================================
# SECTION 1: SAMPLE IDENTIFICATION
# ============================================================================
SAMPLE_NAME <- "<YOUR_SAMPLE_NAME>"  # e.g., "E7.5_rep1"; must match combined pipeline sample name
PROJECT_NAME <- "<YOUR_PROJECT_NAME>_ATAC_only"

# ============================================================================
# SECTION 2: INPUT DATA PATHS (FROM COMBINED PIPELINE)
# ============================================================================
## Path to outputs from the combined pipeline
COMBINED_PIPELINE_DIR <- "<PATH_TO_COMBINED_PIPELINE_OUTPUT>"  # e.g., "/path/to/combined/pipeline/output"

# Seurat object with splits (from combined pipeline Step 02)
INPUT_SEURAT_SPLITS <- file.path(COMBINED_PIPELINE_DIR, "splits", SAMPLE_NAME,
                                  paste0(SAMPLE_NAME, "_seurat_obj_with_splits.rds"))

# Target gene directory (from combined pipeline)
INPUT_TARGET_GENES_DIR <- file.path(COMBINED_PIPELINE_DIR, "target_genes", SAMPLE_NAME)

# ============================================================================
# SECTION 3: OUTPUT DIRECTORIES
# ============================================================================
BASE_OUTPUT_DIR <- "<PATH_TO_OUTPUT_DIRECTORY>"  # e.g., "~/scMultiPreDICT_output/atac_only/processed/"
OUTPUT_METACELLS_DIR <- file.path(BASE_OUTPUT_DIR, "metacells", SAMPLE_NAME)
OUTPUT_FEATURES_DIR <- file.path(BASE_OUTPUT_DIR, "features_extracted", SAMPLE_NAME)

# ============================================================================
# SECTION 4: SPECIES/GENOME CONFIGURATION
# ============================================================================
SPECIES <- "<SPECIES>"  # "mouse" or "human"

if (SPECIES == "mouse") {
  GENOME <- "mm10"
  ENSDB_PACKAGE <- "EnsDb.Mmusculus.v79"
  BSGENOME_PACKAGE <- "BSgenome.Mmusculus.UCSC.mm10"
}
if (SPECIES == "human") {
  GENOME <- "hg38"
  ENSDB_PACKAGE <- "EnsDb.Hsapiens.v86"
  BSGENOME_PACKAGE <- "BSgenome.Hsapiens.UCSC.hg38"
}

# ============================================================================
# SECTION 5: DIMENSIONALITY REDUCTION METHOD
# ============================================================================
# Options: "LSI" (default, fast), "PeakVI" (requires Python)
DIM_REDUCTION_METHOD <- "LSI"

DIMRED_METHOD_SUFFIX <- switch(
  DIM_REDUCTION_METHOD,
  "LSI" = "lsi",
  "PeakVI" = "peakvi",
  tolower(DIM_REDUCTION_METHOD)
)

# Update output dirs with method suffix
OUTPUT_METACELLS_DIR <- file.path(BASE_OUTPUT_DIR, "metacells", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
OUTPUT_FEATURES_DIR <- file.path(BASE_OUTPUT_DIR, "features_extracted", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)

# ============================================================================
# SECTION 6: METACELL/KNN SMOOTHING PARAMETERS
# ============================================================================
LSI_DIMS <- 30
K_NEIGHBORS <- 20
SEED_METACELL <- 2025

# ============================================================================
# SECTION 7: TARGET GENE CONFIGURATION
# ============================================================================
# Use the same target genes as the combined pipeline for fair comparison
HVG_GENE_FILE <- file.path(INPUT_TARGET_GENES_DIR, "<HVG_GENE_FILE>")  # e.g., "target_genes_hvg_100.txt"
RANDOM_GENE_FILE <- file.path(INPUT_TARGET_GENES_DIR, "<RANDOM_GENE_FILE>")  # e.g., "target_genes_random_100.txt"

# Which gene sets to analyze
MODEL_GENE_SET <- "both"  # Options: "HVG", "Random_genes", "both"

# ============================================================================
# SECTION 8: FEATURE EXTRACTION PARAMETERS
# ============================================================================
# For ATAC-only, features are peaks within a genomic window of the target gene
GENE_WINDOW_KB <- 250  # ±250kb window from TSS
MIN_PEAKS_PER_GENE <- 1

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

# Neural network grid/search settings (can be overridden)
NN_USE_GRID_SEARCH <- TRUE
NN_GRID_UNITS <- c(256)
NN_GRID_DROPOUT <- c(0, 0.1)
NN_GRID_BATCH <- c(128, 256)


# Conda environment for TensorFlow/Keras
CONDA_ENV_NAME <- "<YOUR_CONDA_ENV>"

# ============================================================================
# SECTION 10: COMPUTATIONAL RESOURCES & SLURM SETTINGS
# ============================================================================
N_CORES <- 4
MAX_CORES_TRAINING <- 32

CONDA_ENV_R <- "<YOUR_R_CONDA_ENV>"
CONDA_ENV_PYTHON <- "<YOUR_PYTHON_CONDA_ENV>"

SLURM_PARTITION <- "<YOUR_SLURM_PARTITION>"  # e.g., "compute"

SLURM_RESOURCES <- list(
  step_030 = list(cpus = 32, mem_gb = 128, time_hours = 6, job_name = "metacell"),
  step_040 = list(cpus = 32, mem_gb = 128, time_hours = 8, job_name = "feature_extract"),
  step_050 = list(cpus = 32, mem_gb = 200, time_hours = 12, job_name = "model_linear"),
  step_060 = list(cpus = 16, mem_gb = 200, time_hours = 24, job_name = "model_nn")
)

SLURM_LOG_DIR <- "<PATH_TO_SLURM_LOG_DIR>"  # e.g., "~/scATAC_only/jobs/logs"

# ============================================================================
# SECTION 11: RESULTS OUTPUT DIRECTORIES
# ============================================================================
BASE_RESULTS_DIR <- "<PATH_TO_RESULTS_DIRECTORY>"  # e.g., "~/scMultiPreDICT_output/atac_only/results/"

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

print_config <- function() {
  cat("\n=======================================================================\n")
  cat("ATAC-ONLY PIPELINE CONFIGURATION\n")
  cat("=======================================================================\n\n")
  cat("Sample:", SAMPLE_NAME, "\n")
  cat("Species:", SPECIES, "(", GENOME, ")\n")
  cat("Dimensionality Reduction:", DIM_REDUCTION_METHOD, "\n")
  cat("Gene window: ±", GENE_WINDOW_KB, "kb\n")
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

cat("\n[CONFIG] Loaded ATAC-only configuration for sample:", SAMPLE_NAME, "\n")
cat("[CONFIG] Species:", SPECIES, "| Genome:", GENOME, "\n\n")

# ============================================================================
# SBATCH GENERATION (RNA-config style)
# ============================================================================

generate_sbatch <- function(step_name, script_name, output_dir = ".",
                            conda_env = NULL, extra_args = "",
                            array_range = NULL, is_python = FALSE) {

  if (!step_name %in% names(SLURM_RESOURCES)) {
    stop(sprintf("Unknown step: %s. Available: %s",
                 step_name, paste(names(SLURM_RESOURCES), collapse = ", ")))
  }
  res <- SLURM_RESOURCES[[step_name]]

  if (is.null(conda_env)) {
    conda_env <- if (is_python) CONDA_ENV_PYTHON else CONDA_ENV_R
  }

  hours <- floor(res$time_hours)
  minutes <- round((res$time_hours - hours) * 60)
  time_str <- sprintf("%02d:%02d:00", hours, minutes)

  log_dir <- file.path(path.expand(SLURM_LOG_DIR), res$job_name)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  use_gpu <- !is.null(res$gpu) && isTRUE(res$gpu)
  gpu_count <- if (!is.null(res$gpu_count)) res$gpu_count else 1

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

  if (use_gpu) sbatch_lines <- c(sbatch_lines, sprintf("#SBATCH --gres=gpu:%d", gpu_count))
  if (!is.null(array_range)) sbatch_lines <- c(sbatch_lines, sprintf("#SBATCH --array=%s", array_range))

  sbatch_lines <- c(sbatch_lines, "", "set -euo pipefail", "")

  sbatch_lines <- c(sbatch_lines,
    "# Initialize and activate conda environment",
    "source ~/miniconda3/etc/profile.d/conda.sh",
    sprintf("conda activate %s", conda_env),
    ""
  )

  if (!is_python) {
    sbatch_lines <- c(sbatch_lines,
      "# Set R library paths (conda R)",
      "export R_LIBS=\"$CONDA_PREFIX/lib/R/library\"",
      "export R_LIBS_USER=\"$CONDA_PREFIX/lib/R/library\"",
      "export R_LIBS_SITE=\"\"",
      ""
    )
  }

  run_cmd <- if (is_python) {
    sprintf("python %s %s", script_name, extra_args)
  } else {
    sprintf("Rscript %s %s", script_name, extra_args)
  }

  sbatch_lines <- c(sbatch_lines,
    "# Run the script",
    run_cmd,
    "",
    "echo \"Job completed at $(date)\"",
    "echo \"Exit code: $?\""
  )

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  sbatch_filename <- file.path(output_dir, sprintf("%s_%s.sbatch", step_name, SAMPLE_NAME))
  writeLines(sbatch_lines, sbatch_filename)

  cat(sprintf("Generated: %s\n", sbatch_filename))
  cat(sprintf("  Resources: %d CPUs, %dG RAM, %s time\n", res$cpus, res$mem_gb, time_str))
  invisible(sbatch_filename)
}

generate_all_sbatch <- function(output_dir = ".") {
  cat("\n=== Generating SLURM sbatch files (scATAC-only) ===\n\n")

  generate_sbatch("step_030", "03_metacell_creation.R", output_dir)
  generate_sbatch("step_040", "04_feature_extraction.R", output_dir)
  generate_sbatch("step_050", "05_linear_tree_models.R", output_dir)
  generate_sbatch("step_060", "06_neural_network.R", output_dir)
 
  cat("\nAll sbatch files generated in:", output_dir, "\n\n")
  cat("Submission order:\n")
  cat(sprintf("sbatch step_030_%s.sbatch\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_040_%s.sbatch\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_050_%s.sbatch\n", SAMPLE_NAME))
  cat(sprintf("sbatch step_060_%s.sbatch\n", SAMPLE_NAME))
  cat("\n")
}
