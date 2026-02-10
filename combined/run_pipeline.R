#!/usr/bin/env Rscript
# ============================================================================
# scMultiPreDICT - Pipeline Orchestration Script
# ============================================================================
#
# Description:
#   Master script for executing the complete scMultiPreDICT analysis pipeline.
#   Supports selective execution of individual steps or step ranges.
#
# Pipeline Steps:
#   1. Quality Control - Cell filtering based on RNA/ATAC metrics
#   2. Data Splitting - Stratified train/validation/test partitioning (uses pre-computed target gene lists)
#   3. Metacell Creation - k-NN smoothing for noise reduction
#   4. Feature Extraction - Gene-specific feature matrix construction
#   5. Model Training - Linear regression and Random Forest models
#   6. Neural Network Training - Deep learning model training
# 
#
# Usage:
#   Rscript run_pipeline.R                    # Execute all steps
#   Rscript run_pipeline.R --start 3          # Begin from step 3
#   Rscript run_pipeline.R --stop 4           # Stop after step 4
#   Rscript run_pipeline.R --steps 1,2,5,6    # Execute specific steps
#   Rscript run_pipeline.R --help             # Display usage information
#
# ============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default settings
config_path <- "config.R"
start_step <- 1
stop_step <- 6
specific_steps <- NULL

# Argument parsing
i <- 1
while (i <= length(args)) {
  if (args[i] == "--config") {
    config_path <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--start") {
    start_step <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--stop") {
    stop_step <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--steps") {
    specific_steps <- as.integer(strsplit(args[i + 1], ",")[[1]])
    i <- i + 2
  } else if (args[i] == "--help") {
    cat("
scMultiPreDICT Pipeline Orchestration Script

Usage:
  Rscript run_pipeline.R [options]

Options:
  --config PATH    Path to configuration file (default: config.R)
  --start N        Begin execution from step N (default: 1)
  --stop N         Stop execution after step N (default: 6)
  --steps N,M,...  Execute only specified steps
  --help           Display this help message

Pipeline Steps:
  1: Quality Control
  2: Data Splitting
  3: Metacell Creation
  4: Feature Extraction
  5: Linear and Tree-Based Model Training
  6: Neural Network Training

Example:
  Rscript run_pipeline.R --config my_config.R --steps 1,2,3
")
    quit(status = 0)
  } else {
    i <- i + 1
  }
}

# Display header
cat("\n")
cat("============================================================\n")
cat("                    scMultiPreDICT                          \n")
cat("                                                            \n")
cat("   Gene Expression Prediction from Multiome Data            \n")
cat("============================================================\n")
cat("\n")

# Load configuration
if (!file.exists(config_path)) {
  stop("Configuration file not found: ", config_path, 
       "\nPlease copy R/config_template.R to config.R and modify as needed.")
}

source(config_path)
print_config()

# Determine steps to execute
if (!is.null(specific_steps)) {
  steps_to_run <- specific_steps
} else {
  steps_to_run <- seq(start_step, stop_step)
}

cat(sprintf("\nSteps to execute: %s\n", paste(steps_to_run, collapse = ", ")))

# Select metacell creation script based on dimensionality reduction method
metacell_script <- switch(

  DIM_REDUCTION_METHOD,
  "pca_lsi" = "R/03a_metacell_creation_pca_lsi.R",
  "wnn" = "R/03a_metacell_creation_wnn.R",
  "multivi" = "R/03a_metacell_creation_multivi.R",
  "scvi_peakvi" = "R/03a_metacell_creation_scvi_peakvi.R",
  stop(sprintf("Unknown DIM_REDUCTION_METHOD: %s. Valid options: pca_lsi, wnn, multivi, scvi_peakvi", 
               DIM_REDUCTION_METHOD))
)

cat(sprintf("Dimensionality reduction method: %s\n", DIM_REDUCTION_METHOD))
cat(sprintf("Metacell script: %s\n", metacell_script))

# Define pipeline steps
step_info <- list(
  list(name = "Quality Control", 
       script = "R/01_quality_control.R"),
  list(name = "Data Splitting", 
       script = "R/02a_data_splitting.R"),
  list(name = sprintf("Metacell Creation (%s)", DIM_REDUCTION_METHOD), 
       script = metacell_script),
  list(name = "Feature Extraction", 
       script = "R/04_feature_extraction.R"),
  list(name = "Linear and Tree-Based Model Training", 
       script = "R/05_linear_tree_models.R"),
  list(name = "Neural Network Training", 
       script = "R/06_neural_network.R")
)

# Initialize timing
start_time <- Sys.time()
step_times <- list()

# Execute pipeline steps
for (step_num in steps_to_run) {
  
  if (step_num < 1 || step_num > length(step_info)) {
    warning(sprintf("Invalid step number: %d. Skipping.", step_num))
    next
  }
  
  step <- step_info[[step_num]]
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat(sprintf("  STEP %d: %s\n", step_num, step$name))
  cat("------------------------------------------------------------\n")
  
  step_start <- Sys.time()
  
  if (file.exists(step$script)) {
    tryCatch({
      source(step$script, local = new.env())
      step_times[[step_num]] <- difftime(Sys.time(), step_start, units = "mins")
      cat(sprintf("\n[SUCCESS] Step %d completed in %.1f minutes\n", 
                  step_num, as.numeric(step_times[[step_num]])))
    }, error = function(e) {
      cat(sprintf("\n[ERROR] Step %d failed: %s\n", step_num, e$message))
      cat("\nPipeline terminated due to error.\n")
      quit(status = 1)
    })
  } else {
    warning(sprintf("Script not found: %s", step$script))
  }
}

# Display summary
end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("============================================================\n")
cat("                   PIPELINE COMPLETE                        \n")
cat("============================================================\n")
cat("\n")

cat("Execution Time Summary:\n")
for (step_num in names(step_times)) {
  step <- step_info[[as.integer(step_num)]]
  cat(sprintf("  Step %s (%s): %.1f minutes\n", 
              step_num, step$name, as.numeric(step_times[[step_num]])))
}
cat(sprintf("\nTotal execution time: %.1f minutes (%.1f hours)\n", 
            as.numeric(total_time), as.numeric(total_time) / 60))

cat("\nOutput Locations:\n")
cat(sprintf("  Seurat objects:  %s\n", OUTPUT_SEURAT_DIR))
cat(sprintf("  Data splits:     %s\n", OUTPUT_SPLITS_DIR))
cat(sprintf("  Metacells:       %s\n", OUTPUT_METACELLS_DIR))
cat(sprintf("  Features:        %s\n", OUTPUT_FEATURES_DIR))
cat(sprintf("  Model results:   %s\n", OUTPUT_MODELS_LINEAR_DIR))
cat("\n")
