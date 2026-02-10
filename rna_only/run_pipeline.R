#!/usr/bin/env Rscript
# ============================================================================
# scMultiPreDICT - Master Pipeline Runner
# ============================================================================
#
# This script runs the complete scMultiPreDICT pipeline:
#   1. Quality Control
#   2. Data Splitting & Target Gene Selection
#   3. Metacell Creation
#   4. Feature Extraction
#   5. Model Training
#   6. Results Aggregation
#   7. Plotting
#
# Usage:
#   Rscript run_pipeline.R                    # Run all steps
#   Rscript run_pipeline.R --start 3          # Start from step 3
#   Rscript run_pipeline.R --stop 4           # Stop after step 4
#   Rscript run_pipeline.R --steps 1,2,5,6    # Run specific steps
#
# ============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default settings
config_path <- "config.R"
start_step <- 1
stop_step <- 4
specific_steps <- NULL

# Parse arguments
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
scMultiPreDICT Pipeline Runner

Usage:
  Rscript run_pipeline.R [options]

Options:
  --config PATH    Path to config file (default: config.R)
  --start N        Start from step N (default: 1)
  --stop N         Stop after step N (default: 4)
  --steps N,M,...  Run only specific steps
  --help           Show this help message

Steps:
  1: Metacell Creation (PCA)
  2: Feature Extraction
  3: Model Training (Linear + RF)
  4: Neural Network Training

Example:
  Rscript run_pipeline.R --config my_config.R --steps 1,2,3
")
    quit(status = 0)
  } else {
    i <- i + 1
  }
}

# Print header
cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║                                                              ║\n")
cat("║              scMultiPreDICT Pipeline Runner                  ║\n")
cat("║                                                              ║\n")
cat("║   Gene Expression Prediction from Multiome Data             ║\n")
cat("║                                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Load configuration
if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path, 
       "\nPlease copy config_template.R to config.R and edit it.")
}

source(config_path)
print_config()

# Determine which steps to run
if (!is.null(specific_steps)) {
  steps_to_run <- specific_steps
} else {
  steps_to_run <- seq(start_step, stop_step)
}

cat(sprintf("\nSteps to run: %s\n", paste(steps_to_run, collapse = ", ")))

# Define step information
step_info <- list(
  list(name = "Metacell Creation", 
       script = "R/03a_metacell_creation.R"),
  list(name = "Feature Extraction", 
       script = "R/04_feature_extraction.R"),
  list(name = "Model Training", 
       script = "R/05_linear_tree_models.R"),
  list(name = "Neural Network Training", 
       script = "R/06_neural_network.R")
)

# Track timing
start_time <- Sys.time()
step_times <- list()

# Run each step
for (step_num in steps_to_run) {
  
  if (step_num < 1 || step_num > length(step_info)) {
    warning(sprintf("Invalid step number: %d. Skipping.", step_num))
    next
  }
  
  step <- step_info[[step_num]]
  
  cat("\n")
  cat("┌──────────────────────────────────────────────────────────────┐\n")
  cat(sprintf("│  STEP %d: %-52s│\n", step_num, step$name))
  cat("└──────────────────────────────────────────────────────────────┘\n")
  
  step_start <- Sys.time()
  
  # Source the script
  if (file.exists(step$script)) {
    tryCatch({
      source(step$script, local = new.env())
      step_times[[step_num]] <- difftime(Sys.time(), step_start, units = "mins")
      cat(sprintf("\n✓ Step %d completed in %.1f minutes\n", 
                  step_num, as.numeric(step_times[[step_num]])))
    }, error = function(e) {
      cat(sprintf("\n✗ Step %d FAILED: %s\n", step_num, e$message))
      cat("\nPipeline stopped due to error.\n")
      quit(status = 1)
    })
  } else {
    warning(sprintf("Script not found: %s", step$script))
  }
}

# Print summary
end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║                    PIPELINE COMPLETE                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("--- Timing Summary ---\n")
for (step_num in names(step_times)) {
  step <- step_info[[as.integer(step_num)]]
  cat(sprintf("  Step %s (%s): %.1f minutes\n", 
              step_num, step$name, as.numeric(step_times[[step_num]])))
}
cat(sprintf("\nTotal time: %.1f minutes (%.1f hours)\n", 
            as.numeric(total_time), as.numeric(total_time) / 60))

cat("\n--- Output Locations ---\n")
cat(sprintf("  Metacells:   %s\n", OUTPUT_METACELLS_DIR))
cat(sprintf("  Features:    %s\n", OUTPUT_FEATURES_DIR))
cat(sprintf("  Models:      %s\n", OUTPUT_MODELS_LINEAR_DIR))
cat(sprintf("  Results:     %s\n", OUTPUT_MODELS_LINEAR_DIR))

cat("\n")
