
# ============================================================================
# Step_050: Gene Expression Prediction - Linear Models & Random Forest (scRNA-only, Automated)
# ============================================================================
# This script trains multiple ML models to predict gene expression from HVG features only (RNA-only, no peak features).
# Models: OLS, Ridge, Lasso, Elastic Net, Random Forest
# Input: Gene-specific features from Step_040
# Output: Model predictions, metrics, and feature importance
# Usage:
#   1. Ensure config.R is properly configured
#   2. Run: Rscript 05_linear_tree_models.R
# ============================================================================


if (!exists("CONFIG_LOADED")) {
  # Get the directory of this script to source config.R
  get_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      return(dirname(normalizePath(sub("^--file=", "", file_arg))))
    }
    return(".")
  }
  script_dir <- get_script_dir()

  # Source configuration file
  config_path <- file.path(script_dir, "config.R")
  if (!file.exists(config_path)) {
    config_path <- "config.R"
  }
  if (!file.exists(config_path)) {
    stop("config.R not found! Please ensure config.R is in the same directory as this script.")
  }
  cat("Loading configuration from:", config_path, "\n")
  source(config_path)
}


# =========================================================================
# LOAD REQUIRED LIBRARIES
# =========================================================================
suppressPackageStartupMessages({
  library(glmnet)
  library(ranger)
  library(dplyr)
  library(Matrix)
  library(parallel)
  library(doParallel)
  library(ggplot2)
  library(tidyr)
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
cat("STEP 050: Model Training (Linear & Tree-Based, RNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Configuration:\n")
cat(sprintf("  Sample: %s\n", SAMPLE_NAME))
cat(sprintf("  Gene set to process: %s\n", GENE_SET))
cat(sprintf("  Random seed: %d\n", RANDOM_SEED))
cat(sprintf("  Random Forest trees: %d\n", N_TREES))
cat(sprintf("  Parallel cores: %d\n", NUM_CORES))


# ============================================================================
# PUBLICATION-READY PLOTTING SETUP
# ============================================================================

# Colorblind-friendly palette (Wong 2011)
colorblind_palette <- c(
  "#E69F00",  # Orange
  "#56B4E9",  # Sky blue
  "#009E73",  # Bluish green
  "#F0E442",  # Yellow
  "#0072B2",  # Blue
  "#D55E00",  # Vermillion
  "#CC79A7",  # Reddish purple
  "#999999"   # Gray
)

# Model colors
model_colors <- c(
  "OLS" = "#0072B2",           # Blue
  "Ridge" = "#56B4E9",         # Sky blue
  "Lasso" = "#009E73",         # Green
  "ElasticNet" = "#E69F00",    # Orange
  "RandomForest" = "#D55E00"   # Vermillion
)

# Gene set colors
geneset_colors <- c(
  "HVG" = "#0072B2",
  "Random_genes" = "#E69F00"
)

# Publication theme
theme_pub <- function(base_size = 12, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5, 
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 2, hjust = 1, color = "gray50"),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
      strip.text = element_text(size = base_size, face = "bold", color = "black"),
      plot.margin = margin(15, 15, 15, 15)
    )
}

save_plot <- function(plot, filename, width = 10, height = 8, dpi = 600, 
                                  create_pdf = TRUE) {
  ggsave(
    filename = paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  cat(sprintf("  Saved: %s.png\n", basename(filename)))
  
  if (create_pdf) {
    ggsave(
      filename = paste0(filename, ".pdf"),
      plot = plot,
      width = width,
      height = height,
      device = cairo_pdf,
      bg = "white"
    )
    cat(sprintf("  Saved: %s.pdf\n", basename(filename)))
  }
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Min-max scaling
scale_mm <- function(M, mn, rg, zv = integer(0)) {
  M <- as.matrix(M)
  M <- sweep(M, 2, mn, "-")
  M <- sweep(M, 2, rg, "/")
  if (length(zv)) M[, zv] <- 0
  return(M)
}

# Evaluation metrics
rmse <- function(y, p) {
  sqrt(mean((y - p)^2, na.rm = TRUE))
}

r2 <- function(y, p) {
  ss_res <- sum((y - p)^2, na.rm = TRUE)
  ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  return(1 - ss_res / ss_tot)
}

sp <- function(y, p) {
  suppressWarnings(cor(y, p, method = "spearman", use = "complete.obs"))
}

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

train_models_for_gene <- function(gene_data, gene_name, gene_set_name, output_base_dir, 
                                  n_trees = 500, random_seed = 123) {
  
  # Extract data
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
  
  # Scale features
  mins   <- apply(as.matrix(X_tr), 2, min, na.rm = TRUE)
  maxs   <- apply(as.matrix(X_tr), 2, max, na.rm = TRUE)
  ranges <- maxs - mins
  zerov <- which(!is.finite(ranges) | ranges == 0)
  if (length(zerov)) ranges[zerov] <- 1
  
  X_tr <- scale_mm(X_tr, mins, ranges, zerov)
  X_va <- scale_mm(X_va, mins, ranges, zerov)
  X_te <- scale_mm(X_te, mins, ranges, zerov)
  
  rownames(X_tr) <- cells_tr
  rownames(X_va) <- cells_va
  rownames(X_te) <- cells_te
  
  # Train models
  set.seed(random_seed)
  
  # OLS
  df_tr <- data.frame(y = y_tr, as.data.frame(X_tr), check.names = FALSE)
  ols <- tryCatch(lm(y ~ ., data = df_tr), error = function(e) NULL)
  
  # Ridge
  rid <- tryCatch(
    cv.glmnet(X_tr, y_tr, alpha = 0, family = "gaussian", standardize = FALSE),
    error = function(e) NULL
  )
  
  # Lasso
  las <- tryCatch(
    cv.glmnet(X_tr, y_tr, alpha = 1, family = "gaussian", standardize = FALSE),
    error = function(e) NULL
  )
  
  # Elastic Net
  ene <- tryCatch(
    cv.glmnet(X_tr, y_tr, alpha = 0.5, family = "gaussian", standardize = FALSE),
    error = function(e) NULL
  )
  
  # Random Forest
  rf <- tryCatch(
    ranger(
      x = X_tr, y = y_tr,
      mtry = max(1, floor(ncol(X_tr)/3)),
      num.trees = n_trees,
      seed = random_seed,
      importance = "impurity"
    ),
    error = function(e) NULL
  )
  
  # Generate predictions
  make_preds <- function(model_type, model_obj) {
    if (is.null(model_obj)) {
      return(list(
        train = rep(NA_real_, nrow(X_tr)),
        val = rep(NA_real_, nrow(X_va)),
        test = rep(NA_real_, nrow(X_te))
      ))
    }
    
    if (model_type == "OLS") {
      tr <- pred_df(model_obj, X_tr)
      va <- pred_df(model_obj, X_va)
      te <- pred_df(model_obj, X_te)
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
  
  # Calculate metrics
  calc_metrics <- function(name, preds) {
    dplyr::bind_rows(
      data.frame(Model = name, Split = "Train", 
                 RMSE = rmse(y_tr, preds$train), 
                 R2 = r2(y_tr, preds$train), 
                 Spearman = sp(y_tr, preds$train)),
      data.frame(Model = name, Split = "Val", 
                 RMSE = rmse(y_va, preds$val), 
                 R2 = r2(y_va, preds$val), 
                 Spearman = sp(y_va, preds$val)),
      data.frame(Model = name, Split = "Test", 
                 RMSE = rmse(y_te, preds$test), 
                 R2 = r2(y_te, preds$test), 
                 Spearman = sp(y_te, preds$test))
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
  all_metrics$N_Features <- gene_data$n_hvg_features
  all_metrics$Feature_Type <- "RNA_only"
  
  # Save results
  gene_output_dir <- file.path(output_base_dir, gene_name)
  dir.create(gene_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(all_metrics, file.path(gene_output_dir, "metrics.csv"), row.names = FALSE)
  
  # Save predictions
  preds_train <- data.frame(
    cell = cells_tr,
    actual = y_tr,
    OLS = ols_preds$train,
    Ridge = ridge_preds$train,
    Lasso = lasso_preds$train,
    ElasticNet = enet_preds$train,
    RandomForest = rf_preds$train
  )
  write.csv(preds_train, file.path(gene_output_dir, "predictions_train.csv"), row.names = FALSE)
  
  preds_val <- data.frame(
    cell = cells_va,
    actual = y_va,
    OLS = ols_preds$val,
    Ridge = ridge_preds$val,
    Lasso = lasso_preds$val,
    ElasticNet = enet_preds$val,
    RandomForest = rf_preds$val
  )
  write.csv(preds_val, file.path(gene_output_dir, "predictions_val.csv"), row.names = FALSE)
  
  preds_test <- data.frame(
    cell = cells_te,
    actual = y_te,
    OLS = ols_preds$test,
    Ridge = ridge_preds$test,
    Lasso = lasso_preds$test,
    ElasticNet = enet_preds$test,
    RandomForest = rf_preds$test
  )
  write.csv(preds_test, file.path(gene_output_dir, "predictions_test.csv"), row.names = FALSE)
  
  # Save RF importance
  if (!is.null(rf)) {
    rf_importance <- data.frame(
      Feature = colnames(X_tr),
      Importance = rf$variable.importance,
      Original_Feature = original_features
    )
    rf_importance <- rf_importance[order(rf_importance$Importance, decreasing = TRUE), ]
    write.csv(rf_importance, file.path(gene_output_dir, "rf_importance.csv"), row.names = FALSE)
  }
  
  return(all_metrics)
}

# ============================================================================
# PRINT CONFIGURATION
# ============================================================================
print_config()
print_output_directories()

# Create output directories
create_output_directories()

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 070: Model Training (scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat(sprintf("Dimensionality Reduction Method: %s\n", DIM_REDUCTION_METHOD))
cat(sprintf("Output directory: %s\n\n", OUTPUT_MODELS_LINEAR_DIR))

# Parallel processing setup
NUM_CORES <- min(parallel::detectCores() - 1, MAX_CORES_TRAINING)
NUM_CORES <- max(1, NUM_CORES)
cat(sprintf("Using %d cores for parallel processing\n\n", NUM_CORES))

set.seed(SEED_MODEL)

# ============================================================================
# DETERMINE WHICH GENE SETS TO PROCESS
# ============================================================================
cat("=== Determining gene sets to process ===\n\n")

gene_sets_to_process <- c()
if (MODEL_GENE_SET %in% c("HVG", "both")) {
  hvg_features_file <- file.path(OUTPUT_FEATURES_DIR, "HVG", "gene_specific_features.rds")
  if (file.exists(hvg_features_file)) {
    gene_sets_to_process <- c(gene_sets_to_process, "HVG")
    cat("  ✓ HVG gene set found\n")
  } else {
    cat("  ✗ HVG gene set not found:", hvg_features_file, "\n")
  }
}
if (MODEL_GENE_SET %in% c("Random_genes", "both")) {
  random_features_file <- file.path(OUTPUT_FEATURES_DIR, "Random_genes", "gene_specific_features.rds")
  if (file.exists(random_features_file)) {
    gene_sets_to_process <- c(gene_sets_to_process, "Random_genes")
    cat("  ✓ Random_genes gene set found\n")
  } else {
    cat("  ✗ Random_genes gene set not found:", random_features_file, "\n")
  }
}
if (length(gene_sets_to_process) == 0) {
  stop("No gene sets found to process! Run Step_060 first.")
}
cat(sprintf("\nWill process %d gene set(s): %s\n", 
            length(gene_sets_to_process), paste(gene_sets_to_process, collapse = ", ")))

# ============================================================================
# PROCESS EACH GENE SET
# ============================================================================

all_results <- list()

for (gene_set in gene_sets_to_process) {
  
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat(sprintf("PROCESSING GENE SET: %s\n", gene_set))
  cat("=", rep("=", 70), "\n\n", sep = "")
  
  # Load features
  features_file <- file.path(OUTPUT_FEATURES_DIR, gene_set, "gene_specific_features.rds")
  gene_features <- readRDS(features_file)
  
  cat(sprintf("Loaded %d genes from %s\n", length(gene_features), gene_set))
  
  # Create output directory for this gene set
  set_output_dir <- file.path(OUTPUT_MODELS_LINEAR_DIR, gene_set)
  if (!dir.exists(set_output_dir)) {
    dir.create(set_output_dir, recursive = TRUE)
  }
  
  # Train models for all genes in parallel
  cat(sprintf("\nTraining models for %d genes...\n", length(gene_features)))
  
  # Register parallel backend
  cl <- makeCluster(NUM_CORES)
  registerDoParallel(cl)
  
  # Export necessary functions and packages to workers
  clusterExport(cl, c("scale_mm", "rmse", "r2", "sp", "pred_df", "train_models_for_gene",
                      "RF_N_TREES", "SEED_MODEL", "gene_features"))
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(glmnet)
      library(ranger)
      library(dplyr)
    })
  })
  
  # Train models in parallel
  gene_names <- names(gene_features)
  
  # Prepare a list of arguments for each gene
  args_list <- lapply(gene_names, function(gene_name) {
    list(
      gene_name = gene_name,
      gene_data = gene_features[[gene_name]],
      gene_set = gene_set,
      set_output_dir = set_output_dir,
      RF_N_TREES = RF_N_TREES,
      SEED_MODEL = SEED_MODEL
    )
  })
  
  results_list <- parLapply(cl, args_list, function(args) {
    tryCatch({
      train_models_for_gene(
        gene_data = args$gene_data,
        gene_name = args$gene_name,
        gene_set_name = args$gene_set,
        output_base_dir = args$set_output_dir,
        n_trees = args$RF_N_TREES,
        random_seed = args$SEED_MODEL
      )
    }, error = function(e) {
      return(data.frame(
        Gene = args$gene_name,
        Error = as.character(e$message),
        stringsAsFactors = FALSE
      ))
    })
  })
  
  stopCluster(cl)
  
  # Combine results
  results_df <- dplyr::bind_rows(results_list)
  
  # Filter out error rows
  if ("Error" %in% colnames(results_df)) {
    errors <- results_df[!is.na(results_df$Error), ]
    if (nrow(errors) > 0) {
      cat(sprintf("\nWarning: %d genes failed to train:\n", nrow(errors)))
      print(head(errors))
    }
    results_df <- results_df[is.na(results_df$Error) | !("Error" %in% colnames(results_df)), ]
  }
  
  # Save aggregated results
  agg_file <- file.path(set_output_dir, "aggregated_metrics.csv")
  write.csv(results_df, agg_file, row.names = FALSE)
  cat(sprintf("\nSaved aggregated metrics: %s\n", agg_file))
  
  all_results[[gene_set]] <- results_df
  
  cat(sprintf("\n✓ Completed %s gene set: %d genes processed\n", gene_set, length(gene_names)))
}

# ============================================================================
# GENERATE SUMMARY PLOTS
# ============================================================================
cat("\n=== Generating summary plots ===\n\n")

plots_dir <- file.path(OUTPUT_MODELS_LINEAR_DIR, "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Combine all results
combined_results <- dplyr::bind_rows(all_results, .id = "GeneSet_ID")

if (nrow(combined_results) > 0 && "R2" %in% colnames(combined_results)) {
  
  # Filter to test set only for main plots
  test_results <- combined_results %>%
    filter(Split == "Test")
  
  # ----------------------------------------------------------------------------
  # 1. R² DISTRIBUTION BY MODEL (Test Set)
  # ----------------------------------------------------------------------------
  p_r2_box <- ggplot(test_results, aes(x = Model, y = R2, fill = Model)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = model_colors) +
    labs(
      x = "Model",
      y = expression(R^2),
      title = "Model Performance Comparison (Test Set)",
      subtitle = sprintf("%s | %s | n=%d genes", SAMPLE_NAME, DIMRED_METHOD_SUFFIX, 
                         length(unique(test_results$Gene)))
    ) +
    theme_pub() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  save_plot(p_r2_box, 
            file.path(plots_dir, paste0(SAMPLE_NAME, "_r2_by_model_boxplot")),
                        width = 10, height = 6)
  
  # ----------------------------------------------------------------------------
  # 2. R² DISTRIBUTION BY MODEL AND GENE SET
  # ----------------------------------------------------------------------------
  if (length(unique(test_results$GeneSet)) > 1) {
    p_r2_facet <- ggplot(test_results, aes(x = Model, y = R2, fill = Model)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      scale_fill_manual(values = model_colors) +
      facet_wrap(~ GeneSet) +
      labs(
        x = "Model",
        y = expression(R^2),
        title = "Model Performance by Gene Set (Test Set)",
        subtitle = sprintf("%s | %s", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
      ) +
      theme_pub() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    save_plot(p_r2_facet, 
              file.path(plots_dir, paste0(SAMPLE_NAME, "_r2_by_model_geneset")),
              width = 14, height = 6)
  }
  
  # ----------------------------------------------------------------------------
  # 3. SPEARMAN CORRELATION BY MODEL
  # ----------------------------------------------------------------------------
  p_sp_box <- ggplot(test_results, aes(x = Model, y = Spearman, fill = Model)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = model_colors) +
    labs(
      x = "Model",
      y = "Spearman Correlation",
      title = "Spearman Correlation by Model (Test Set)",
      subtitle = sprintf("%s | %s", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
    ) +
    theme_pub() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  save_plot(p_sp_box, 
            file.path(plots_dir, paste0(SAMPLE_NAME, "_spearman_by_model_boxplot")),
            width = 10, height = 6)
  
  # ----------------------------------------------------------------------------
  # 4. TRAIN VS TEST R² (Overfitting Check)
  # ----------------------------------------------------------------------------
  train_test_comparison <- combined_results %>%
    filter(Split %in% c("Train", "Test")) %>%
    select(Gene, GeneSet, Model, Split, R2) %>%
    pivot_wider(names_from = Split, values_from = R2, names_prefix = "R2_")
  
  p_overfit <- ggplot(train_test_comparison, aes(x = R2_Train, y = R2_Test, color = Model)) +
    geom_point(alpha = 0.4, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = model_colors) +
    labs(
      x = expression(R^2 ~ "(Train)"),
      y = expression(R^2 ~ "(Test)"),
      title = "Train vs Test Performance (Overfitting Check)",
      subtitle = sprintf("%s | Points below diagonal = overfitting", SAMPLE_NAME)
    ) +
    theme_pub() +
    coord_fixed()
  
  save_plot(p_overfit, 
            file.path(plots_dir, paste0(SAMPLE_NAME, "_train_vs_test_r2")),
            width = 10, height = 10)
  
  # ----------------------------------------------------------------------------
  # 5. BEST MODEL PER GENE
  # ----------------------------------------------------------------------------
  best_per_gene <- test_results %>%
    group_by(Gene, GeneSet) %>%
    slice_max(order_by = R2, n = 1) %>%
    ungroup()
  
  best_model_counts <- best_per_gene %>%
    count(Model, GeneSet) %>%
    mutate(Model = factor(Model, levels = names(model_colors)))
  
  p_best_model <- ggplot(best_model_counts, aes(x = Model, y = n, fill = Model)) +
    geom_col(color = "black", alpha = 0.8) +
    geom_text(aes(label = n), vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_manual(values = model_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      x = "Model",
      y = "Number of Genes",
      title = "Best Model per Gene (by Test R²)",
      subtitle = sprintf("%s | %s", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
    ) +
    theme_pub() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  save_plot(p_best_model, 
            file.path(plots_dir, paste0(SAMPLE_NAME, "_best_model_per_gene")),
            width = 10, height = 6)
  
  # ----------------------------------------------------------------------------
  # 6. R² HISTOGRAM
  # ----------------------------------------------------------------------------
  p_r2_hist <- ggplot(test_results, aes(x = R2, fill = Model)) +
    geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = model_colors) +
    facet_wrap(~ Model, ncol = 5) +
    labs(
      x = expression(R^2),
      y = "Number of Genes",
      title = "Distribution of Test R² by Model",
      subtitle = sprintf("%s | %s", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
    ) +
    theme_pub() +
    theme(legend.position = "none")
  
  save_plot(p_r2_hist, 
            file.path(plots_dir, paste0(SAMPLE_NAME, "_r2_histogram_by_model")),
            width = 14, height = 5)
  
  # ----------------------------------------------------------------------------
  # 7. SUMMARY STATISTICS TABLE
  # ----------------------------------------------------------------------------
  summary_stats <- test_results %>%
    group_by(Model) %>%
    summarise(
      Mean_R2 = mean(R2, na.rm = TRUE),
      Median_R2 = median(R2, na.rm = TRUE),
      SD_R2 = sd(R2, na.rm = TRUE),
      Mean_Spearman = mean(Spearman, na.rm = TRUE),
      Mean_RMSE = mean(RMSE, na.rm = TRUE),
      N_Positive_R2 = sum(R2 > 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary_file <- file.path(plots_dir, paste0(SAMPLE_NAME, "_model_summary_stats.csv"))
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat(sprintf("Saved summary statistics: %s\n", summary_file))
  
  # Print summary
  cat("\n=== MODEL PERFORMANCE SUMMARY (Test Set) ===\n\n")
  print(summary_stats)
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 050 COMPLETE (scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Model training completed for all gene sets\n\n")

cat("Models trained:\n")
cat("  - OLS (Ordinary Least Squares)\n")
cat("  - Ridge Regression (L2)\n")
cat("  - Lasso Regression (L1)\n")
cat("  - Elastic Net (L1 + L2)\n")
cat("  - Random Forest\n")

cat("\nOutput files per gene:\n")
cat("  - metrics.csv (train/val/test performance)\n")
cat("  - predictions_train.csv\n")
cat("  - predictions_val.csv\n")
cat("  - predictions_test.csv\n")
cat("  - rf_importance.csv (Random Forest feature importance)\n")

cat("\nAggregated results:\n")
for (gene_set in gene_sets_to_process) {
  cat(sprintf("  - %s/aggregated_metrics.csv\n", gene_set))
}

cat("\nplots saved in plots/ directory\n")

cat("\nNext: Run 06_neural_network.R\n")
