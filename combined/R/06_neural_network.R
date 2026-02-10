# ============================================================================
# Step_060: Gene Expression Prediction - Deep Neural Network (Automated)
# ============================================================================
# This script trains deep neural networks to predict gene expression from
# ATAC peaks + HVG features. It reads configuration from config.R.
#
# Model Architecture:
#   - Configurable hidden layers (1-3)
#   - ReLU activation
#   - Optional dropout regularization
#   - Adam optimizer with early stopping
#
# Input: Gene-specific features from Step_040
# Output: Model predictions, metrics, and trained Keras models
#
# Requirements:
#   - Conda environment with TensorFlow/Keras installed
#   - Configured in config.R via CONDA_ENV_NAME
#
# Usage:
#   1. Ensure config.R is properly configured
#   2. Run: Rscript 06_neural_network.R
# ============================================================================

# Record start time
start_time <- Sys.time()

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
# PYTHON / TENSORFLOW ENVIRONMENT
# ============================================================================
suppressMessages(library(reticulate))

# Use configured conda environment
use_condaenv(CONDA_ENV_NAME, required = TRUE)
cat("Python configuration:\n")
print(py_config())

# Import TensorFlow and Keras
tf <- reticulate::import("tensorflow", delay_load = TRUE)
keras <- tf$keras

# CPU-friendly thread settings for cluster computing
Sys.setenv(TF_NUM_INTRAOP_THREADS = "1", TF_NUM_INTEROP_THREADS = "1")

# ============================================================================
# LOAD REMAINING LIBRARIES
# ============================================================================
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(matrixStats)
  library(readr)
})

# Create output directories
create_output_directories()

# Set parameters from config
GENE_SET <- MODEL_GENE_SET
RANDOM_SEED <- SEED_MODEL

set.seed(RANDOM_SEED)

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 060: Model Training (Deep Neural Network)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Configuration:\n")
cat(sprintf("  Sample: %s\n", SAMPLE_NAME))
cat(sprintf("  Gene set to process: %s\n", GENE_SET))
cat(sprintf("  Random seed: %d\n", RANDOM_SEED))
cat(sprintf("  Hidden layers: %d\n", NN_N_HIDDEN_LAYERS))
cat(sprintf("  Hidden units: %d\n", NN_HIDDEN_UNITS))
cat(sprintf("  Dropout rate: %.2f\n", NN_DROPOUT_RATE))
cat(sprintf("  Learning rate: %.4f\n", NN_LEARNING_RATE))
cat(sprintf("  Use grid search: %s\n", NN_USE_GRID_SEARCH))

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Min-max scaling for features
scale_mm <- function(M, mn, rg, zv = integer(0)) {
  M <- as.matrix(M)
  M <- sweep(M, 2, mn, "-")
  M <- sweep(M, 2, rg, "/")
  if (length(zv)) M[, zv] <- 0
  return(M)
}

# Root Mean Squared Error
rmse <- function(y, p) {
  sqrt(mean((y - p)^2, na.rm = TRUE))
}

# R-squared
r2 <- function(y, p) {
  ss_res <- sum((y - p)^2, na.rm = TRUE)
  ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  return(1 - ss_res / ss_tot)
}

# Spearman correlation
sp <- function(y, p) {
  suppressWarnings(cor(y, p, method = "spearman", use = "complete.obs"))
}

# Build Keras model with configurable hidden layers
build_keras_model <- function(input_dim, n_units = 256, n_hidden_layers = 3,
                              dropout_rate = 0.1, learning_rate = 1e-3) {
  
  model <- tf$keras$models$Sequential()
  
  # First hidden layer (with input shape)
  model$add(
    tf$keras$layers$Dense(
      units = as.integer(n_units),
      activation = "relu",
      input_shape = list(as.integer(input_dim))
    )
  )
  if (dropout_rate > 0) {
    model$add(tf$keras$layers$Dropout(rate = dropout_rate))
  }
  
  # Second hidden layer (if n_hidden_layers >= 2)
  if (n_hidden_layers >= 2) {
    model$add(
      tf$keras$layers$Dense(
        units = as.integer(n_units / 2),
        activation = "relu"
      )
    )
    if (dropout_rate > 0) {
      model$add(tf$keras$layers$Dropout(rate = dropout_rate))
    }
  }
  
  # Third hidden layer (if n_hidden_layers >= 3)
  if (n_hidden_layers >= 3) {
    model$add(
      tf$keras$layers$Dense(
        units = as.integer(n_units / 4),
        activation = "relu"
      )
    )
    if (dropout_rate > 0) {
      model$add(tf$keras$layers$Dropout(rate = dropout_rate))
    }
  }
  
  # Output layer
  model$add(
    tf$keras$layers$Dense(
      units = as.integer(1),
      activation = "linear"
    )
  )
  
  # Compile with Adam + MSE
  model$compile(
    optimizer = tf$keras$optimizers$Adam(learning_rate = learning_rate),
    loss = "mse",
    metrics = list("mse")
  )
  
  return(model)
}

# ============================================================================
# MAIN TRAINING FUNCTION FOR ONE GENE
# ============================================================================

train_dl_for_gene <- function(gene_data, gene_name, gene_set_name, output_base_dir) {
  
  cat(sprintf("\n=== Processing Gene: %s ===\n", gene_name))
  
  # -------------------------
  # STEP 1: EXTRACT DATA
  # -------------------------
  X_tr <- gene_data$train$X
  y_tr <- gene_data$train$y
  cells_tr <- gene_data$train$cells
  
  X_va <- gene_data$validation$X
  y_va <- gene_data$validation$y
  cells_va <- gene_data$validation$cells
  
  X_te <- gene_data$test$X
  y_te <- gene_data$test$y
  cells_te <- gene_data$test$cells
  
  original_features <- colnames(X_tr)
  
  n_total_features <- ncol(X_tr)
  n_peaks <- ifelse(!is.null(gene_data$n_peaks), gene_data$n_peaks, 0)
  n_hvg_features <- ifelse(!is.null(gene_data$n_hvg_features), gene_data$n_hvg_features, 0)
  
  cat(sprintf("  Features: %d total\n", n_total_features))
  cat(sprintf("  Train: %d cells, Val: %d cells, Test: %d cells\n",
              nrow(X_tr), nrow(X_va), nrow(X_te)))
  
  # -------------------------
  # STEP 1.5: REMOVE TARGET GENE FROM PREDICTORS
  # -------------------------
  hvg_feature_names <- colnames(X_tr)
  target_in_features <- grepl(paste0("^", gene_name, "$"), hvg_feature_names)
  
  if (any(target_in_features)) {
    cat(sprintf("  ⚠ Target gene '%s' found in predictors - removing it\n", gene_name))
    
    X_tr <- X_tr[, !target_in_features, drop = FALSE]
    X_va <- X_va[, !target_in_features, drop = FALSE]
    X_te <- X_te[, !target_in_features, drop = FALSE]
    
    original_features <- original_features[!target_in_features]
    
    n_removed <- sum(target_in_features)
    cat(sprintf("    Removed %d feature(s). New total: %d features\n", n_removed, ncol(X_tr)))
  }
  
  # Skip if no features remain
  if (ncol(X_tr) == 0) {
    cat("  ✗ No predictors left after filtering. Skipping gene.\n")
    return(list(gene = gene_name, status = "failed", error = "No predictors after filtering"))
  }
  
  # -------------------------
  # STEP 2: SCALE FEATURES
  # -------------------------
  cat("  → Scaling features...\n")
  
  mins <- apply(as.matrix(X_tr), 2, min, na.rm = TRUE)
  maxs <- apply(as.matrix(X_tr), 2, max, na.rm = TRUE)
  ranges <- maxs - mins
  
  zerov <- which(!is.finite(ranges) | ranges == 0)
  
  if (length(zerov)) {
    ranges[zerov] <- 1
    cat(sprintf("    Warning: %d zero-variance features detected\n", length(zerov)))
  }
  
  X_tr <- scale_mm(X_tr, mins, ranges, zerov)
  X_va <- scale_mm(X_va, mins, ranges, zerov)
  X_te <- scale_mm(X_te, mins, ranges, zerov)
  
  rownames(X_tr) <- cells_tr
  rownames(X_va) <- cells_va
  rownames(X_te) <- cells_te
  
  # -------------------------
  # STEP 3: PREPARE TENSORS
  # -------------------------
  x_tr <- as.matrix(X_tr)
  x_va <- as.matrix(X_va)
  x_te <- as.matrix(X_te)
  
  y_tr_vec <- as.numeric(y_tr)
  y_va_vec <- as.numeric(y_va)
  y_te_vec <- as.numeric(y_te)
  
  # Check for zero variance in labels
  if (sd(y_tr_vec, na.rm = TRUE) == 0) {
    cat("  ✗ Zero variance in training labels. Skipping gene.\n")
    return(list(gene = gene_name, status = "failed", error = "Zero variance in training labels"))
  }
  
  # Standardize target variable for training stability
  y_mu <- mean(y_tr_vec, na.rm = TRUE)
  y_sd <- sd(y_tr_vec, na.rm = TRUE)
  if (!is.finite(y_sd) || y_sd == 0) y_sd <- 1
  
  y_tr_s <- (y_tr_vec - y_mu) / y_sd
  y_va_s <- (y_va_vec - y_mu) / y_sd
  
  # Convert to TensorFlow tensors
  x_tr_tf <- tf$convert_to_tensor(x_tr, dtype = tf$float32)
  y_tr_tf <- tf$convert_to_tensor(matrix(y_tr_s, ncol = 1), dtype = tf$float32)
  
  x_va_tf <- tf$convert_to_tensor(x_va, dtype = tf$float32)
  y_va_tf <- tf$convert_to_tensor(matrix(y_va_s, ncol = 1), dtype = tf$float32)
  
  x_te_tf <- tf$convert_to_tensor(x_te, dtype = tf$float32)
  
  # -------------------------
  # STEP 4: TRAIN MODEL (with or without grid search)
  # -------------------------
  
  if (NN_USE_GRID_SEARCH) {
    cat("  → Tuning deep learning model (grid search)...\n")
    
    grid <- expand.grid(
      n_units = NN_GRID_UNITS,
      dropout_rate = NN_GRID_DROPOUT,
      learning_rate = c(NN_LEARNING_RATE),
      batch_size = NN_GRID_BATCH,
      stringsAsFactors = FALSE
    )
    
    results_grid <- list()
    best_val_RMSE <- Inf
    best_model <- NULL
    best_params <- NULL
    
    for (i in seq_len(nrow(grid))) {
      params <- grid[i, ]
      
      cat(sprintf("    - Combo %d/%d: units=%d, dropout=%.2f, batch=%d\n",
                  i, nrow(grid), params$n_units, params$dropout_rate, params$batch_size))
      
      model <- build_keras_model(
        input_dim = ncol(x_tr),
        n_units = params$n_units,
        n_hidden_layers = NN_N_HIDDEN_LAYERS,
        dropout_rate = params$dropout_rate,
        learning_rate = params$learning_rate
      )
      
      history <- tryCatch({
        model$fit(
          x = x_tr_tf,
          y = y_tr_tf,
          validation_data = reticulate::tuple(x_va_tf, y_va_tf),
          batch_size = as.integer(params$batch_size),
          epochs = as.integer(NN_MAX_EPOCHS),
          verbose = 0L,
          callbacks = list(
            tf$keras$callbacks$EarlyStopping(
              monitor = "val_loss",
              patience = as.integer(NN_EARLY_STOP_PATIENCE),
              restore_best_weights = TRUE
            )
          )
        )
      }, error = function(e) {
        cat("      ✗ Training failed: ", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(history)) {
        results_grid[[i]] <- cbind(params, val_RMSE = NA_real_, val_R2 = NA_real_)
        next
      }
      
      # Validation predictions
      va_pred_s <- as.numeric(model$predict(x_va_tf, batch_size = as.integer(1024)))
      va_pred <- va_pred_s * y_sd + y_mu
      
      val_RMSE <- rmse(y_va_vec, va_pred)
      val_R2 <- r2(y_va_vec, va_pred)
      
      results_grid[[i]] <- cbind(params, val_RMSE = val_RMSE, val_R2 = val_R2)
      
      if (!is.na(val_RMSE) && val_RMSE < best_val_RMSE) {
        best_val_RMSE <- val_RMSE
        best_model <- model
        best_params <- params
      }
    }
    
    # Use last model if none succeeded
    if (is.null(best_model)) {
      best_model <- model
      best_params <- params
      best_val_RMSE <- ifelse(is.na(val_RMSE), Inf, val_RMSE)
    }
    
  } else {
    # Fixed parameters (no grid search)
    cat("  → Training deep learning model (fixed parameters)...\n")
    
    best_params <- data.frame(
      n_units = NN_HIDDEN_UNITS,
      dropout_rate = NN_DROPOUT_RATE,
      learning_rate = NN_LEARNING_RATE,
      batch_size = NN_BATCH_SIZE
    )
    
    best_model <- build_keras_model(
      input_dim = ncol(x_tr),
      n_units = NN_HIDDEN_UNITS,
      n_hidden_layers = NN_N_HIDDEN_LAYERS,
      dropout_rate = NN_DROPOUT_RATE,
      learning_rate = NN_LEARNING_RATE
    )
    
    history <- tryCatch({
      best_model$fit(
        x = x_tr_tf,
        y = y_tr_tf,
        validation_data = reticulate::tuple(x_va_tf, y_va_tf),
        batch_size = as.integer(NN_BATCH_SIZE),
        epochs = as.integer(NN_MAX_EPOCHS),
        verbose = 0L,
        callbacks = list(
          tf$keras$callbacks$EarlyStopping(
            monitor = "val_loss",
            patience = as.integer(NN_EARLY_STOP_PATIENCE),
            restore_best_weights = TRUE
          )
        )
      )
    }, error = function(e) {
      cat("  ✗ Training failed: ", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(history)) {
      return(list(gene = gene_name, status = "failed", error = "Training failed"))
    }
    
    results_grid <- list()
    va_pred_s <- as.numeric(best_model$predict(x_va_tf, batch_size = as.integer(1024)))
    va_pred <- va_pred_s * y_sd + y_mu
    best_val_RMSE <- rmse(y_va_vec, va_pred)
  }
  
  cat("  → Best DL params:\n")
  print(best_params)
  
  # -------------------------
  # STEP 5: PREDICTIONS & METRICS
  # -------------------------
  cat("  → Generating predictions and metrics...\n")
  
  pred_tr_s <- as.numeric(best_model$predict(x_tr_tf, batch_size = as.integer(1024)))
  pred_va_s <- as.numeric(best_model$predict(x_va_tf, batch_size = as.integer(1024)))
  pred_te_s <- as.numeric(best_model$predict(x_te_tf, batch_size = as.integer(1024)))
  
  # Rescale to original scale
  pred_tr <- pred_tr_s * y_sd + y_mu
  pred_va <- pred_va_s * y_sd + y_mu
  pred_te <- pred_te_s * y_sd + y_mu
  
  all_metrics <- dplyr::bind_rows(
    data.frame(Model = "DeepNN", Split = "Train", 
               RMSE = rmse(y_tr_vec, pred_tr), R2 = r2(y_tr_vec, pred_tr), 
               Spearman = sp(y_tr_vec, pred_tr)),
    data.frame(Model = "DeepNN", Split = "Val", 
               RMSE = rmse(y_va_vec, pred_va), R2 = r2(y_va_vec, pred_va), 
               Spearman = sp(y_va_vec, pred_va)),
    data.frame(Model = "DeepNN", Split = "Test", 
               RMSE = rmse(y_te_vec, pred_te), R2 = r2(y_te_vec, pred_te), 
               Spearman = sp(y_te_vec, pred_te))
  )
  
  all_metrics$Gene <- gene_name
  all_metrics$GeneSet <- gene_set_name
  all_metrics$N_Features <- ncol(x_tr)
  all_metrics$N_Peaks <- n_peaks
  all_metrics$N_HVG <- n_hvg_features
  
  # -------------------------
  # STEP 6: SAVE OUTPUTS
  # -------------------------
  gene_output_dir <- file.path(output_base_dir, gene_name)
  dir.create(gene_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("  → Saving results...\n")
  
  # Save best model
  model_path <- normalizePath(file.path(gene_output_dir, "best_model.keras"), mustWork = FALSE)
  best_model$save(model_path)
  
  # Save scaled feature matrices
  saveRDS(list(X_train = X_tr, X_val = X_va, X_test = X_te), 
          file.path(gene_output_dir, "scaled_features.rds"))
  
  # Metrics
  write_csv(all_metrics, file.path(gene_output_dir, "metrics.csv"))
  
  # Predictions
  write_csv(data.frame(Cell_ID = cells_tr, Actual = y_tr_vec, DeepNN = pred_tr),
            file.path(gene_output_dir, "predictions_train.csv"))
  write_csv(data.frame(Cell_ID = cells_va, Actual = y_va_vec, DeepNN = pred_va),
            file.path(gene_output_dir, "predictions_val.csv"))
  write_csv(data.frame(Cell_ID = cells_te, Actual = y_te_vec, DeepNN = pred_te),
            file.path(gene_output_dir, "predictions_test.csv"))
  
  # Grid search results
  if (NN_USE_GRID_SEARCH && length(results_grid) > 0) {
    grid_df <- do.call(rbind, results_grid)
    write_csv(grid_df, file.path(gene_output_dir, "grid_search_results.csv"))
  }
  
  # Metadata
  metadata <- data.frame(
    Gene = gene_name,
    GeneSet = gene_set_name,
    N_Total_Features = n_total_features,
    N_Peaks = n_peaks,
    N_HVG_Features = n_hvg_features,
    N_Train_Cells = nrow(x_tr),
    N_Val_Cells = nrow(x_va),
    N_Test_Cells = nrow(x_te),
    N_Zero_Variance_Features = length(zerov),
    N_Hidden_Layers = NN_N_HIDDEN_LAYERS,
    Best_n_units = best_params$n_units,
    Best_dropout_rate = best_params$dropout_rate,
    Best_batch_size = best_params$batch_size,
    Best_val_RMSE = best_val_RMSE,
    RandomSeed = RANDOM_SEED,
    Analysis_Date = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  write_csv(metadata, file.path(gene_output_dir, "analysis_metadata.csv"))
  
  cat("  ✓ Saved DL outputs for gene:", gene_name, "\n")
  
  return(list(gene = gene_name, status = "success", metrics = all_metrics, metadata = metadata))
}

# ============================================================================
# AGGREGATE RESULTS FUNCTION
# ============================================================================

aggregate_results <- function(results, output_dir, gene_set_name) {
  
  cat(sprintf("\n=== Aggregating results for %s ===\n", gene_set_name))
  
  if (length(results) == 0) {
    cat("No results to aggregate.\n")
    return(NULL)
  }
  
  metrics <- dplyr::bind_rows(lapply(results, function(x) {
    if (x$status == "success") x$metrics else NULL
  }))
  
  if (nrow(metrics) > 0) {
    write_csv(metrics, file.path(output_dir, "all_genes_combined_metrics.csv"))
    cat(sprintf("  ✓ Saved combined metrics for %d genes\n", dplyr::n_distinct(metrics$Gene)))
    
    summary <- metrics %>%
      group_by(Model, Split) %>%
      summarise(
        N_Genes = n_distinct(Gene),
        Mean_R2 = mean(R2, na.rm = TRUE),
        Median_R2 = median(R2, na.rm = TRUE),
        SD_R2 = sd(R2, na.rm = TRUE),
        Mean_RMSE = mean(RMSE, na.rm = TRUE),
        Mean_Spearman = mean(Spearman, na.rm = TRUE),
        .groups = "drop"
      )
    
    write_csv(summary, file.path(output_dir, "summary_statistics.csv"))
    cat("  ✓ Saved summary statistics\n")
    
    cat("\nPerformance Summary (Test Set):\n")
    print(summary %>% filter(Split == "Test"), n = Inf)
  }
  
  statuses <- sapply(results, function(x) x$status)
  genes <- sapply(results, function(x) x$gene)
  
  status_report <- data.frame(Gene = genes, Status = statuses, GeneSet = gene_set_name)
  write_csv(status_report, file.path(output_dir, "processing_status.csv"))
  cat("  ✓ Saved processing status report\n")
  
  n_success <- sum(statuses == "success")
  n_failed <- sum(statuses != "success")
  cat(sprintf("\nProcessing Summary: %d succeeded, %d failed\n", n_success, n_failed))
  
  return(status_report)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Helper function to convert number to word for folder naming
number_to_word <- function(n) {
  words <- c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten")
  if (n >= 1 && n <= 10) {
    return(words[n])
  } else {
    return(as.character(n))
  }
}

# Create architecture folder name based on hidden layers
arch_folder <- paste0(number_to_word(NN_N_HIDDEN_LAYERS), "_hidden_layer")
cat(sprintf("\nArchitecture: %s (%d hidden layers)\n", arch_folder, NN_N_HIDDEN_LAYERS))

# Input and output directories
input_dir <- path.expand(OUTPUT_FEATURES_DIR)
output_dir <- file.path(path.expand(OUTPUT_MODELS_NN_DIR), arch_folder)

cat(sprintf("\nInput directory: %s\n", input_dir))
cat(sprintf("Output directory: %s\n", output_dir))

# Auto-detect gene sets by looking for subdirectories with feature files
available_gene_sets <- c()
for (potential_set in c("HVG", "Random_genes")) {
  feature_file <- file.path(input_dir, potential_set, "gene_specific_features.rds")
  if (file.exists(feature_file)) {
    available_gene_sets <- c(available_gene_sets, potential_set)
  }
}

# Filter based on config setting
if (GENE_SET == "HVG") {
  gene_sets_to_process <- intersect(available_gene_sets, "HVG")
} else if (GENE_SET == "Random_genes") {
  gene_sets_to_process <- intersect(available_gene_sets, "Random_genes")
} else {
  # "both" or any other value - process all available
  gene_sets_to_process <- available_gene_sets
}

if (length(gene_sets_to_process) == 0) {
  stop("No gene sets found! Please run Step_060 first to extract features.")
}

cat(sprintf("\nAvailable gene sets: %s\n", paste(available_gene_sets, collapse = ", ")))
cat(sprintf("Gene sets to process: %s\n", paste(gene_sets_to_process, collapse = ", ")))

# Process each gene set
for (gene_set_name in gene_sets_to_process) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf(">>> Processing Gene Set: %s <<<\n", gene_set_name))
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  feature_file <- file.path(input_dir, gene_set_name, "gene_specific_features.rds")
  set_output_dir <- file.path(output_dir, gene_set_name)
  
  if (!file.exists(feature_file)) {
    cat(sprintf("ERROR: Feature file not found: %s\n", feature_file))
    cat("Skipping this gene set.\n")
    next
  }
  
  if (!dir.exists(set_output_dir)) {
    dir.create(set_output_dir, recursive = TRUE)
  }
  
  cat(sprintf("\nLoading features from: %s\n", feature_file))
  gene_features <- readRDS(feature_file)
  cat(sprintf("Loaded features for %d genes\n", length(gene_features)))
  
  results <- list()
  n_genes <- length(gene_features)
  
  cat(sprintf("\nProcessing %d genes...\n", n_genes))
  
  idx <- 0
  for (gene_name in names(gene_features)) {
    idx <- idx + 1
    gene_data <- gene_features[[gene_name]]
    
    cat(sprintf("\n[%d/%d] Gene: %s\n", idx, n_genes, gene_name))
    
    result <- tryCatch({
      train_dl_for_gene(gene_data, gene_name, gene_set_name, set_output_dir)
    }, error = function(e) {
      cat(sprintf("  ✗ ERROR: %s\n", e$message))
      return(list(gene = gene_name, status = "failed", error = e$message))
    })
    
    results[[gene_name]] <- result
  }
  
  cat("\n", paste(rep("-", 70), collapse = ""), "\n")
  status_report <- aggregate_results(results, set_output_dir, gene_set_name)
  
  cat(sprintf("\n>>> Completed %s gene set <<<\n", gene_set_name))
  cat(paste(rep("=", 70), collapse = ""), "\n")
}

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 071 COMPLETE\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat(sprintf("Total runtime: %.2f minutes\n\n", as.numeric(runtime)))

cat("Model trained:\n")
cat(sprintf("  - DeepNN (%d hidden layers, ReLU activation)\n", NN_N_HIDDEN_LAYERS))
cat(sprintf("  - Architecture: %d → %d → %d → 1\n", 
            NN_HIDDEN_UNITS, NN_HIDDEN_UNITS/2, NN_HIDDEN_UNITS/4))
cat(sprintf("  - Dropout: %.1f%%\n", NN_DROPOUT_RATE * 100))

cat("\nMetrics calculated:\n")
cat("  - RMSE (Root Mean Squared Error)\n")
cat("  - R² (Coefficient of Determination)\n")
cat("  - Spearman Correlation\n")

cat("\nOutput directory:", output_dir, "\n\n")

cat("Output structure per gene:\n")
cat("  ├── best_model.keras (trained Keras model)\n")
cat("  ├── metrics.csv\n")
cat("  ├── predictions_train.csv\n")
cat("  ├── predictions_val.csv\n")
cat("  ├── predictions_test.csv\n")
cat("  ├── grid_search_results.csv (if grid search enabled)\n")
cat("  ├── scaled_features.rds\n")
cat("  └── analysis_metadata.csv\n")

cat("\n", paste(rep("=", 70), collapse = ""), "\n\n")
