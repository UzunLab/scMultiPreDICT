# ============================================================================
# Step_050: Gene Expression Prediction - Linear Models & Random Forest (Automated)
# ============================================================================
# This script trains multiple ML models to predict gene expression from
# ATAC peaks + HVG features. It reads configuration from config.R.
#
# Models trained:
#   1. OLS (Ordinary Least Squares)
#   2. Ridge Regression (L2 regularization)
#   3. Lasso Regression (L1 regularization)
#   4. Elastic Net (L1 + L2 regularization)
#   5. Random Forest (ensemble of decision trees)
#
# Input: Gene-specific features from Step_040 
# Output: Model predictions, metrics, and feature importance
#
# Usage:
#   1. Ensure config_template.R is properly configured
#   2. Run: Rscript 05_linear_tree_models.R
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
  library(glmnet)
  library(ranger)
  library(dplyr)
  library(Matrix)
  library(parallel)
  library(doParallel)
})

# Create output directories
create_output_directories()

# Set parameters from config
GENE_SET <- MODEL_GENE_SET
RANDOM_SEED <- SEED_MODEL
N_TREES <- RF_N_TREES
NUM_CORES <- min(parallel::detectCores() - 1, MAX_CORES_TRAINING)
NUM_CORES <- max(1, NUM_CORES)

set.seed(RANDOM_SEED)

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 050: Model Training (Linear & Tree-Based)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Configuration:\n")
cat(sprintf("  Sample: %s\n", SAMPLE_NAME))
cat(sprintf("  Gene set to process: %s\n", GENE_SET))
cat(sprintf("  Random seed: %d\n", RANDOM_SEED))
cat(sprintf("  Random Forest trees: %d\n", N_TREES))
cat(sprintf("  Parallel cores: %d\n", NUM_CORES))

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Scale matrix columns to [0, 1] range using min-max normalization
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

# Generate predictions from a fitted model
pred_df <- function(fit, X) {
  if (is.null(fit)) {
    return(rep(NA_real_, nrow(X)))
  } else {
    return(as.numeric(predict(fit, as.data.frame(X))))
  }
}

# ============================================================================
# MAIN TRAINING FUNCTION FOR ONE GENE
# ============================================================================

train_models_for_gene <- function(gene_data, gene_name, gene_set_name, output_base_dir) {
  
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
  
  cat(sprintf("  Features: %d peaks + %d HVG features = %d total\n", 
              gene_data$n_peaks, gene_data$n_hvg_features, gene_data$n_total_features))
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
  
  # -------------------------
  # STEP 2: SCALE FEATURES
  # -------------------------
  cat("  → Scaling features...\n")
  
  mins   <- apply(as.matrix(X_tr), 2, min, na.rm = TRUE)
  maxs   <- apply(as.matrix(X_tr), 2, max, na.rm = TRUE)
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
  # STEP 3: TRAIN MODELS
  # -------------------------
  cat("  → Training models...\n")
  set.seed(RANDOM_SEED)
  
  # MODEL 1: OLS
  cat("    - OLS... ")
  df_tr <- data.frame(y = y_tr, as.data.frame(X_tr), check.names = FALSE)
  ols <- tryCatch(lm(y ~ ., data = df_tr), error = function(e) { cat("FAILED\n"); NULL })
  if (!is.null(ols)) cat("✓\n")
  
  # MODEL 2: Ridge
  cat("    - Ridge... ")
  rid <- tryCatch(cv.glmnet(X_tr, y_tr, alpha = 0, family = "gaussian", standardize = FALSE),
                  error = function(e) { cat("FAILED\n"); NULL })
  if (!is.null(rid)) cat("✓\n")
  
  # MODEL 3: Lasso
  cat("    - Lasso... ")
  las <- tryCatch(cv.glmnet(X_tr, y_tr, alpha = 1, family = "gaussian", standardize = FALSE),
                  error = function(e) { cat("FAILED\n"); NULL })
  if (!is.null(las)) cat("✓\n")
  
  # MODEL 4: Elastic Net
  cat("    - Elastic Net... ")
  ene <- tryCatch(cv.glmnet(X_tr, y_tr, alpha = 0.5, family = "gaussian", standardize = FALSE),
                  error = function(e) { cat("FAILED\n"); NULL })
  if (!is.null(ene)) cat("✓\n")
  
  # MODEL 5: Random Forest
  cat("    - Random Forest... ")
  rf <- tryCatch(
    ranger(x = X_tr, y = y_tr, mtry = max(1, floor(ncol(X_tr)/3)),
           num.trees = N_TREES, seed = RANDOM_SEED, importance = "impurity"),
    error = function(e) { cat("FAILED\n"); NULL })
  if (!is.null(rf)) cat("✓\n")
  
  # -------------------------
  # STEP 4: GENERATE PREDICTIONS
  # -------------------------
  cat("  → Generating predictions...\n")
  
  make_preds <- function(model_type, model_obj) {
    if (is.null(model_obj)) {
      return(list(train = rep(NA_real_, nrow(X_tr)), val = rep(NA_real_, nrow(X_va)), 
                  test = rep(NA_real_, nrow(X_te))))
    }
    
    if (model_type == "OLS") {
      tr <- pred_df(model_obj, X_tr); va <- pred_df(model_obj, X_va); te <- pred_df(model_obj, X_te)
    } else if (model_type == "RF") {
      tr <- as.numeric(predict(model_obj, data = as.data.frame(X_tr))$predictions)
      va <- as.numeric(predict(model_obj, data = as.data.frame(X_va))$predictions)
      te <- as.numeric(predict(model_obj, data = as.data.frame(X_te))$predictions)
    } else {
      tr <- as.numeric(predict(model_obj, X_tr, s = "lambda.min"))
      va <- as.numeric(predict(model_obj, X_va, s = "lambda.min"))
      te <- as.numeric(predict(model_obj, X_te, s = "lambda.min"))
    }
    return(list(train = tr, val = va, test = te))
  }
  
  ols_preds <- make_preds("OLS", ols)
  ridge_preds <- make_preds("Ridge", rid)
  lasso_preds <- make_preds("Lasso", las)
  enet_preds <- make_preds("ENet", ene)
  rf_preds <- make_preds("RF", rf)
  
  # -------------------------
  # STEP 5: CALCULATE METRICS
  # -------------------------
  cat("  → Calculating metrics...\n")
  
  calc_metrics <- function(name, preds) {
    dplyr::bind_rows(
      data.frame(Model = name, Split = "Train", RMSE = rmse(y_tr, preds$train), 
                 R2 = r2(y_tr, preds$train), Spearman = sp(y_tr, preds$train)),
      data.frame(Model = name, Split = "Val", RMSE = rmse(y_va, preds$val), 
                 R2 = r2(y_va, preds$val), Spearman = sp(y_va, preds$val)),
      data.frame(Model = name, Split = "Test", RMSE = rmse(y_te, preds$test), 
                 R2 = r2(y_te, preds$test), Spearman = sp(y_te, preds$test))
    )
  }
  
  all_metrics <- dplyr::bind_rows(
    calc_metrics("OLS", ols_preds),
    calc_metrics("Ridge", ridge_preds),
    calc_metrics("Lasso", lasso_preds),
    calc_metrics("ElasticNet", enet_preds),
    calc_metrics("RandomForest", rf_preds)
  )
  
  all_metrics$Gene <- gene_name
  all_metrics$GeneSet <- gene_set_name
  all_metrics$N_Features <- gene_data$n_total_features
  all_metrics$N_Peaks <- gene_data$n_peaks
  all_metrics$N_HVG <- gene_data$n_hvg_features
  
  # -------------------------
  # STEP 6: EXTRACT COEFFICIENTS
  # -------------------------
  cat("  → Extracting coefficients...\n")
  
  # OLS coefficients
  if (!is.null(ols)) {
    ols_coef <- data.frame(Feature = names(coef(ols))[-1], Coefficient = as.numeric(coef(ols)[-1]), Model = "OLS")
    ols_coef$Original_Feature <- original_features[match(ols_coef$Feature, colnames(X_tr))]
  } else {
    ols_coef <- data.frame(Feature = character(0), Coefficient = numeric(0), Model = character(0), Original_Feature = character(0))
  }
  
  # Ridge coefficients
  if (!is.null(rid)) {
    ridge_coef <- data.frame(Feature = colnames(X_tr), Coefficient = as.numeric(coef(rid, s = "lambda.min")[-1]), Model = "Ridge")
    ridge_coef$Original_Feature <- original_features
  } else {
    ridge_coef <- data.frame(Feature = character(0), Coefficient = numeric(0), Model = character(0), Original_Feature = character(0))
  }
  
  # Lasso coefficients
  if (!is.null(las)) {
    lasso_coef <- data.frame(Feature = colnames(X_tr), Coefficient = as.numeric(coef(las, s = "lambda.min")[-1]), Model = "Lasso")
    lasso_coef$Original_Feature <- original_features
  } else {
    lasso_coef <- data.frame(Feature = character(0), Coefficient = numeric(0), Model = character(0), Original_Feature = character(0))
  }
  
  # ElasticNet coefficients
  if (!is.null(ene)) {
    enet_coef <- data.frame(Feature = colnames(X_tr), Coefficient = as.numeric(coef(ene, s = "lambda.min")[-1]), Model = "ElasticNet")
    enet_coef$Original_Feature <- original_features
  } else {
    enet_coef <- data.frame(Feature = character(0), Coefficient = numeric(0), Model = character(0), Original_Feature = character(0))
  }
  
  # Random Forest importance
  if (!is.null(rf)) {
    rf_importance <- data.frame(Feature = colnames(X_tr), Importance = rf$variable.importance, Model = "RandomForest")
    rf_importance$Original_Feature <- original_features
    rf_importance <- rf_importance[order(rf_importance$Importance, decreasing = TRUE), ]
  } else {
    rf_importance <- data.frame(Feature = character(0), Importance = numeric(0), Model = character(0), Original_Feature = character(0))
  }
  
  # Top 20 features
  all_coefs <- rbind(
    if (nrow(ols_coef) > 0) ols_coef[, c("Original_Feature", "Coefficient", "Model")] else NULL,
    if (nrow(ridge_coef) > 0) ridge_coef[, c("Original_Feature", "Coefficient", "Model")] else NULL,
    if (nrow(lasso_coef) > 0) lasso_coef[, c("Original_Feature", "Coefficient", "Model")] else NULL,
    if (nrow(enet_coef) > 0) enet_coef[, c("Original_Feature", "Coefficient", "Model")] else NULL
  )
  
  if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
    all_coefs$Abs_Coefficient <- abs(all_coefs$Coefficient)
    top20_per_model <- all_coefs %>%
      group_by(Model) %>%
      arrange(desc(Abs_Coefficient)) %>%
      slice_head(n = 20) %>%
      ungroup()
  } else {
    top20_per_model <- data.frame(Original_Feature = character(0), Coefficient = numeric(0), 
                                  Model = character(0), Abs_Coefficient = numeric(0))
  }
  
  # -------------------------
  # STEP 7: SAVE OUTPUTS
  # -------------------------
  gene_output_dir <- file.path(output_base_dir, gene_name)
  dir.create(gene_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("  → Saving results...\n")
  
  write.csv(all_metrics, file.path(gene_output_dir, "metrics.csv"), row.names = FALSE)
  
  # Predictions
  pred_train <- data.frame(Cell_ID = cells_tr, Actual = y_tr, OLS = ols_preds$train,
                           Ridge = ridge_preds$train, Lasso = lasso_preds$train,
                           ElasticNet = enet_preds$train, RandomForest = rf_preds$train)
  write.csv(pred_train, file.path(gene_output_dir, "predictions_train.csv"), row.names = FALSE)
  
  pred_val <- data.frame(Cell_ID = cells_va, Actual = y_va, OLS = ols_preds$val,
                         Ridge = ridge_preds$val, Lasso = lasso_preds$val,
                         ElasticNet = enet_preds$val, RandomForest = rf_preds$val)
  write.csv(pred_val, file.path(gene_output_dir, "predictions_val.csv"), row.names = FALSE)
  
  pred_test <- data.frame(Cell_ID = cells_te, Actual = y_te, OLS = ols_preds$test,
                          Ridge = ridge_preds$test, Lasso = lasso_preds$test,
                          ElasticNet = enet_preds$test, RandomForest = rf_preds$test)
  write.csv(pred_test, file.path(gene_output_dir, "predictions_test.csv"), row.names = FALSE)
  
  # Coefficients
  write.csv(ols_coef, file.path(gene_output_dir, "coefficients_ols.csv"), row.names = FALSE)
  write.csv(ridge_coef, file.path(gene_output_dir, "coefficients_ridge.csv"), row.names = FALSE)
  write.csv(lasso_coef, file.path(gene_output_dir, "coefficients_lasso.csv"), row.names = FALSE)
  write.csv(enet_coef, file.path(gene_output_dir, "coefficients_enet.csv"), row.names = FALSE)
  write.csv(top20_per_model, file.path(gene_output_dir, "coefficients_top20.csv"), row.names = FALSE)
  write.csv(rf_importance, file.path(gene_output_dir, "rf_importance.csv"), row.names = FALSE)
  
  # Trained models
  if (!is.null(ols)) saveRDS(ols, file.path(gene_output_dir, "model_ols.rds"))
  if (!is.null(rid)) saveRDS(rid, file.path(gene_output_dir, "model_ridge.rds"))
  if (!is.null(las)) saveRDS(las, file.path(gene_output_dir, "model_lasso.rds"))
  if (!is.null(ene)) saveRDS(ene, file.path(gene_output_dir, "model_elasticnet.rds"))
  if (!is.null(rf)) saveRDS(rf, file.path(gene_output_dir, "model_rf.rds"))
  
  # Scaling parameters (needed to preprocess new data for prediction)
  saveRDS(list(mins = mins, ranges = ranges, zerov = zerov), 
          file.path(gene_output_dir, "scaling_params.rds"))
          
  # Metadata
  metadata <- data.frame(
    Gene = gene_name, GeneSet = gene_set_name,
    N_Total_Features = gene_data$n_total_features,
    N_Peaks = gene_data$n_peaks, N_HVG_Features = gene_data$n_hvg_features,
    N_Train_Cells = nrow(X_tr), N_Val_Cells = nrow(X_va), N_Test_Cells = nrow(X_te),
    N_Zero_Variance_Features = length(zerov),
    OLS_Success = !is.null(ols), Ridge_Success = !is.null(rid),
    Lasso_Success = !is.null(las), ElasticNet_Success = !is.null(ene),
    RandomForest_Success = !is.null(rf),
    RandomSeed = RANDOM_SEED, Analysis_Date = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  write.csv(metadata, file.path(gene_output_dir, "analysis_metadata.csv"), row.names = FALSE)
  
  cat(sprintf("  ✓ Saved outputs to: %s\n", gene_output_dir))
  
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
    write.csv(metrics, file.path(output_dir, "all_genes_combined_metrics.csv"), row.names = FALSE)
    cat(sprintf("  ✓ Saved combined metrics for %d genes\n", n_distinct(metrics$Gene)))
    
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
    
    write.csv(summary, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
    cat("  ✓ Saved summary statistics\n")
    
    cat("\nPerformance Summary (Test Set):\n")
    print(summary %>% filter(Split == "Test"), n = Inf)
  }
  
  statuses <- sapply(results, function(x) x$status)
  genes <- sapply(results, function(x) x$gene)
  
  status_report <- data.frame(Gene = genes, Status = statuses, GeneSet = gene_set_name)
  write.csv(status_report, file.path(output_dir, "processing_status.csv"), row.names = FALSE)
  cat("  ✓ Saved processing status report\n")
  
  n_success <- sum(statuses == "success")
  n_failed <- sum(statuses != "success")
  cat(sprintf("\nProcessing Summary: %d succeeded, %d failed\n", n_success, n_failed))
  
  return(status_report)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Input and output directories
input_dir <- path.expand(OUTPUT_FEATURES_DIR)
output_dir <- path.expand(OUTPUT_MODELS_LINEAR_DIR)

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
# Accept either "Random" or "Random_genes" as synonyms for the random gene set
if (GENE_SET == "HVG") {
  gene_sets_to_process <- intersect(available_gene_sets, "HVG")
} else if (GENE_SET == "Random_genes") {
  gene_sets_to_process <- intersect(available_gene_sets, "Random_genes")
} else {
  # "both" or any other value - process all available
  gene_sets_to_process <- available_gene_sets
}

if (length(gene_sets_to_process) == 0) {
  stop("No gene sets found! Please run Step_040 first to extract features.")
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
  
  n_genes <- length(gene_features)
  cat(sprintf("\nProcessing %d genes using %d cores...\n", n_genes, NUM_CORES))
  
  # Set up parallel cluster
  cat(sprintf("Setting up parallel cluster with %d workers...\n", NUM_CORES))
  cl <- makeCluster(NUM_CORES)
  
  clusterExport(cl, c(
    "train_models_for_gene", "scale_mm", "rmse", "r2", "sp", "pred_df",
    "set_output_dir", "gene_set_name", "RANDOM_SEED", "N_TREES"
  ), envir = environment())
  
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(glmnet)
      library(ranger)
      library(dplyr)
      library(Matrix)
    })
  })
  
  process_one_gene <- function(gene_info) {
    gene_name <- gene_info$name
    gene_data <- gene_info$data
    
    tryCatch(
      train_models_for_gene(gene_data, gene_name, gene_set_name, set_output_dir),
      error = function(e) list(gene = gene_name, status = "failed", error = e$message)
    )
  }
  
  clusterExport(cl, "process_one_gene", envir = environment())
  
  gene_list <- lapply(seq_along(gene_features), function(i) {
    list(name = names(gene_features)[i], data = gene_features[[i]])
  })
  
  cat("Processing genes in parallel...\n")
  start_time <- Sys.time()
  
  results <- parLapply(cl, gene_list, process_one_gene)
  
  stopCluster(cl)
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat(sprintf("\nCompleted processing %d genes in %.1f seconds (%.2f genes/sec)\n", 
              n_genes, elapsed_time, n_genes/elapsed_time))
  
  names(results) <- names(gene_features)
  
  cat("\n", paste(rep("-", 70), collapse = ""), "\n")
  status_report <- aggregate_results(results, set_output_dir, gene_set_name)
  
  cat(sprintf("\n>>> Completed %s gene set <<<\n", gene_set_name))
  cat(paste(rep("=", 70), collapse = ""), "\n")
}

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 070 COMPLETE\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Model training completed successfully!\n\n")

cat("Models trained:\n")
cat("  1. OLS (Ordinary Least Squares)\n")
cat("  2. Ridge Regression (L2 regularization)\n")
cat("  3. Lasso Regression (L1 regularization)\n")
cat("  4. Elastic Net (L1 + L2 regularization)\n")
cat("  5. Random Forest (ensemble of decision trees)\n")

cat("\nMetrics calculated:\n")
cat("  - RMSE (Root Mean Squared Error)\n")
cat("  - R² (Coefficient of Determination)\n")
cat("  - Spearman Correlation\n")

cat("\nOutput directory:", output_dir, "\n\n")

cat("Next step: Run 06_neural_networks.R for deep learning models\n")
