
#============================================================
# Gene-Specific Prediction Evaluation and Feature Importance
# ============================================================
# This script:
# 1) Iterates over multiple datasets and multiple target genes
# 2) Loads gene-specific prediction and feature-importance outputs
# 3) Applies quantile normalization to observed vs predicted expression
# 4) Generates:
# - Observed vs Predicted scatter plots (Spearman correlation)
# - Top-N feature importance bar plots
# 5) Saves outputs in dataset/gene-specific folders
# Designed for:
# - Single-cell multiome gene expression prediction
# - Random Forest models using RNA + ATAC features
# ============================================================

#!/usr/bin/env Rscript

# ==============================================================================
# 1) SETUP: libraries, configuration, global plotting theme
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(preprocessCore)
})


# ------------------------------------------------------------------------------
# Dataset configuration
# Each row corresponds to one dataset.
# - Dataset : short human-readable identifier
# - BasePath : path containing gene-specific model outputs
# - Genes : vector of genes to process for this dataset
# (NA or empty -> auto-select best gene via metrics file)
# ------------------------------------------------------------------------------
# Update BASE_RESULTS_DIR to point to your scMultiPreDICT output directory

BASE_RESULTS_DIR <- "~/scMultiPreDICT_output/results/models"
DIMRED_METHOD <- "pca_lsi"      
GENE_SET <- "HVG"                # Options: "HVG", "Random_genes"

cfg <- tibble::tribble(
  ~Dataset,    ~BasePath,                                                                      ~Genes,
  "E7.5_REP1", file.path(BASE_RESULTS_DIR, "LINEAR_AND_TREE_BASED/E7.5_rep1", DIMRED_METHOD, GENE_SET), c("Etv6", "Stox2", "Mob3b", "Meis2", "Runx1t1"),
  "E7.5_REP2", file.path(BASE_RESULTS_DIR, "LINEAR_AND_TREE_BASED/E7.5_rep2", DIMRED_METHOD, GENE_SET), c("Tbx3", "Bahcc1", "Cask", "Pbx3", "Cbx7"),
  "T_Cells",   file.path(BASE_RESULTS_DIR, "LINEAR_AND_TREE_BASED/T_Cells", DIMRED_METHOD, GENE_SET),   c("RUNX3", "LEF1", "TCF7", "JUN", "STAT4")
)

# Root output directory
out_dir <- "~/scMultiPreDICT_output/paper_figures/gene_specific_figures/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Global Theme (set once) ----
# Ensures visual consistency across all figures
theme_nature <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# ==============================================================================
# 2) HELPERS
# ==============================================================================

# Create sub directories for each gene
make_gene_dir <- function(out_dir, ds, gene) {
  d <- file.path(out_dir, ds, gene)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# Safe CSV reader 
safe_read_csv <- function(path) {
  readr::read_csv(path, show_col_types = FALSE)
}

# Utility: find the first file matching a pattern
# Used for flexible file discovery across model outputs
find_first_file <- function(dir, pattern, recursive = TRUE) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = recursive)
  if (length(f) == 0) return(NA_character_)
  f[1]
}

# Auto-select the best gene based on model performance metrics
# - Prefers Random Forest R² if available
# - Used when no gene list is provided
find_best_gene <- function(base_path) {
  # Prefer a top-level metrics file; if not found, search recursively
  f <- find_first_file(base_path, pattern = "metrics.*\\.csv", recursive = FALSE)
  if (is.na(f)) f <- find_first_file(base_path, pattern = "metrics.*\\.csv", recursive = TRUE)
  if (is.na(f)) return(NA_character_)
  
  metrics <- safe_read_csv(f)
  
  # Normalize target gene column name if needed
  metrics <- metrics %>%
    rename(Target_Gene = any_of(c("Target_Gene", "Gene", "Refgene")))
  
  if (!"Target_Gene" %in% names(metrics)) return(NA_character_)
  if (!"R2" %in% names(metrics)) return(NA_character_)
  
  if ("Model" %in% names(metrics)) {
    best <- metrics %>%
      filter(str_detect(Model, regex("RandomForest|Random Forest|RF", ignore_case = TRUE))) %>%
      arrange(desc(R2)) %>%
      slice(1)
  } else {
    best <- metrics %>% arrange(desc(R2)) %>% slice(1)
  }
  
  if (nrow(best) == 0) return(NA_character_)
  best$Target_Gene[[1]]
}

# Load gene-specific prediction and importance files
# Expects:
# <BasePath>/<Gene>/predictions_test.csv
# <BasePath>/<Gene>/rf_importance.csv
load_gene_files <- function(base_path, gene) {
  gene_dir <- file.path(base_path, gene)
  if (!dir.exists(gene_dir)) stop("Gene folder not found: ", gene_dir)
  
  pred_file <- file.path(gene_dir, "predictions_test.csv")
  if (!file.exists(pred_file)) pred_file <- find_first_file(gene_dir, "predictions.*\\.csv", recursive = FALSE)
  if (is.na(pred_file) || !file.exists(pred_file)) stop("No predictions file found in: ", gene_dir)
  
  imp_file <- file.path(gene_dir, "rf_importance.csv")
  if (!file.exists(imp_file)) imp_file <- find_first_file(gene_dir, "importance.*\\.csv", recursive = FALSE)
  if (is.na(imp_file) || !file.exists(imp_file)) stop("No importance file found in: ", gene_dir)
  
  list(
    preds = data.table::fread(pred_file),
    imps  = data.table::fread(imp_file)
  )
}

# Quantile-normalize two numeric columns
# Implemented without multithreading for HPC compatibility
qn_two_cols <- function(df, x, y) { 
  mat <- as.matrix(df[, c(x, y), drop = FALSE]) 
  
  # Target distribution = mean across sorted columns 
  target <- rowMeans(apply(mat, 2, sort, na.last = TRUE)) 
  qn_col <- function(v) { 
    o <- order(v, na.last = TRUE) 
    out <- v 
    # assign target to sorted positions 
    out[o] <- target 
    out 
  } 
  
  df[[paste0(x, "_QN")]] <- qn_col(mat[, 1]) 
  df[[paste0(y, "_QN")]] <- qn_col(mat[, 2]) 
  df 
  
} 

# Scatter plot: Observed vs Predicted (QN)
# Performance metric: Spearman rank correlation (ρ)
make_scatter_qn <- function(plot_df, gene, ds) {
  
  obs  <- plot_df$Observed_QN
  pred <- plot_df$Predicted_QN
  
  # Spearman rank correlation
  spearman_rho <- suppressWarnings(
    cor(obs, pred, method = "spearman", use = "complete.obs")
  )
  
  ggplot(plot_df, aes(x = Observed_QN, y = Predicted_QN)) +
    geom_point(alpha = 0.35, size = 0.9, color = "#377EB8") +
    geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", linewidth = 0.8, color = "red"
    ) +
    coord_equal() +
    labs(
      title = paste0(gene, " (Quantile Normalized)"),
      subtitle = paste0(ds, " | Random Forest"),
      x = "Observed Expression (QN)",
      y = "Predicted Expression (QN)"
    ) +
    annotate(
      "label",
      x = Inf, y = -Inf,
      label = sprintf("Spearman ρ = %.2f", spearman_rho),
      hjust = 1.05, vjust = -0.2,
      label.size = 0.25,
      size = 4
    ) +
    theme_nature() +
    theme(panel.grid = element_blank())
}

standardize_importance_cols <- function(imps) {
  nm <- tolower(names(imps))
  names(imps) <- nm
  
  feat_col  <- nm[grep("feature|variable|var", nm)][1]
  score_col <- nm[grep("importance|score|incmse|incnode|gain", nm)][1]
  
  if (is.na(feat_col) || is.na(score_col)) return(NULL)
  
  data.table::setnames(imps, c(feat_col, score_col), c("Feature", "Score"))
  imps
}

# Feature type classifier
# RNA features = gene symbols
# ATAC features = genomic coordinates
feature_type <- function(x) {
  # Heuristic: ATAC peaks usually look like chr#:start-end or contain ':'
  ifelse(str_detect(x, regex("^chr|:", ignore_case = TRUE)), "ATAC", "RNA")
}


# Top-N feature importance bar plot with fixed color mapping
make_top20_importance <- function(imps_std, gene, ds) {
  
  top_feats <- as_tibble(imps_std) %>%
    mutate(Score = as.numeric(Score)) %>%
    filter(!is.na(Score)) %>%
    arrange(desc(Score)) %>%
    slice_head(n = 20) %>%
    mutate(
      Type = factor(feature_type(Feature), levels = c("RNA", "ATAC")),
      Feature = fct_reorder(Feature, Score)
    )
  
  ggplot(top_feats, aes(x = Score, y = Feature, fill = Type)) +
    geom_col(width = 0.75, color = "black", linewidth = 0.25) +
    scale_fill_manual(values = c(RNA = "#377EB8", ATAC = "#E41A1C"), drop = FALSE) +
    labs(
      title = paste0(gene, " Top 20 Features"),
      subtitle = paste0(ds, " | Random Forest Importance"),
      x = "Importance Score",
      y = NULL,
      fill = NULL
    ) +
    theme_nature() +
    theme(axis.text.y = element_text(size = 20))
}

# Save plot helper (PDF, consistent resolution)
save_plot <- function(p, path, w, h) {
  ggsave(path, p, width = w, height = h, device = cairo_pdf)
}

# ==============================================================================
# 3) MAIN EXECUTION
# ==============================================================================

purrr::pwalk(cfg, function(Dataset, BasePath, Genes) {
  
  message(sprintf("\nProcessing %s...", Dataset))
  
  # If Genes contains NA or is empty, fall back to best gene
  if (length(Genes) == 0 || all(is.na(Genes)) ) {
    best <- find_best_gene(BasePath)
    if (is.na(best) || best == "") {
      warning(sprintf("  Skipping %s: no gene found from metrics.", Dataset))
      return(invisible(NULL))
    }
    Genes <- best
  }
  
  for (gene in Genes) {
    
    # If a specific entry is NA, replace it with best gene
    if (is.na(gene) || gene == "") {
      gene <- find_best_gene(BasePath)
      if (is.na(gene) || gene == "") {
        warning(sprintf("  %s: could not auto-pick a gene.", Dataset))
        next
      }
    }
    
    message(sprintf("  -> Gene: %s", gene))
    
    gene_dir <- make_gene_dir(out_dir, Dataset, gene)
    
    dat <- tryCatch(
      load_gene_files(BasePath, gene),
      error = function(e) {
        warning(sprintf("  %s | %s: %s", Dataset, gene, e$message))
        return(NULL)
      }
    )
    if (is.null(dat)) next
    
    preds <- dat$preds
    if (!all(c("Actual", "RandomForest") %in% names(preds))) {
      warning(sprintf("  %s | %s: missing Actual/RandomForest in predictions.", Dataset, gene))
      next
    }
    
    plot_df <- preds[, .(Observed = Actual, Predicted = RandomForest)] %>% as.data.frame()
    plot_df <- qn_two_cols(plot_df, "Observed", "Predicted")
    
    # Scatter
    p1 <- make_scatter_qn(plot_df, gene, Dataset)
    save_plot(
      p1,
      file.path(gene_dir, paste0(gene, "_Scatter_QN.pdf")),
      w = 5.2, h = 5.0
    )
    
    # Importance
    imps_std <- standardize_importance_cols(dat$imps)
    if (is.null(imps_std)) {
      warning(sprintf("  %s | %s: could not identify Feature/Score columns.", Dataset, gene))
      next
    }
    
    p2 <- make_top20_importance(imps_std, gene, Dataset)
    save_plot(
      p2,
      file.path(gene_dir, paste0(gene, "_Top20.pdf")),
      w = 8.2, h = 7.0
    )
    
    message(sprintf("  -> Saved: %s/%s", Dataset, gene))
  }
})

message("\n Done! Check folder: ", out_dir)
