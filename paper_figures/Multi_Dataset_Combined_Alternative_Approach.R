#!/usr/bin/env Rscript
# =============================================================================
# Multi-Dataset Combined Performance Plots (Linear/RF + DeepNN)
# =============================================================================
# Publication-oriented plotting of model performance across multiple datasets.
#
# Assumptions:
# - Each dataset has two separate result folders:
#  (1) linear_path: OLS/Ridge/Lasso/ElasticNet/RandomForest
#  (2) nn_path:  DeepNN
# - Each folder contains a combined metrics CSV (common candidates supported)
#
# Outputs:
# - High-res PNG (600 dpi) + vector PDF for each figure
# - Summary CSV with per-dataset, per-model statistics
#
# Usage:
#  Rscript Multi_Dataset_Combined_Alternative_Approach.R
# =============================================================================

start_time <- Sys.time()

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(scales)
  library(Cairo)
})

# =============================================================================
# CONFIGURATION (EDIT)
# =============================================================================
# Update BASE_RESULTS_DIR to point to your scMultiPreDICT output directory

BASE_RESULTS_DIR <- "~/scMultiPreDICT_output/rna_only/results/models"
DIMRED_METHOD <- "pca"       # Options: "pca"
GENE_SET <- "HVG"            # Options: "HVG", "Random_genes"

# Define samples to analyze
SAMPLES <- c("E7.5_rep1", "E7.5_rep2", "T_Cells")
SAMPLE_LABELS <- c("E7.5_REP1", "E7.5_REP2", "T_Cells")

# Build dataset paths
datasets <- setNames(lapply(SAMPLES, function(s) {
  list(
    linear_path = file.path(BASE_RESULTS_DIR, "LINEAR_AND_TREE_BASED", s, DIMRED_METHOD, GENE_SET),
    nn_path = file.path(BASE_RESULTS_DIR, "NEURAL_NETWORKS", s, DIMRED_METHOD, "Three_hidden_layer", GENE_SET)
  )
}), SAMPLE_LABELS)

# Output directory
output_dir <- "~/scMultiPreDICT_output/paper_figures/alternative_approach/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Colorblind-friendly dataset palette (Wong)
dataset_palette <- c(
  "E7.5_REP1" = "#E69F00",
  "E7.5_REP2" = "#56B4E9",
  "T_Cells" = "#009E73"
)

model_levels <- c("OLS", "Ridge", "Lasso", "ElasticNet", "RandomForest", "DeepNN")

# =============================================================================
# PLOTTING THEME (PUBLICATION)
# =============================================================================
theme_publication <- function(base_size = 12, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle  = element_text(size = base_size, hjust = 0.5, color = "gray35"),
      plot.caption = element_text(size = base_size - 2, color = "gray40", hjust = 1),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text  = element_text(size = base_size - 1, color = "black"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text  = element_text(size = base_size - 1),
      panel.grid.major= element_line(color = "gray92", linewidth = 0.3),
      panel.grid.minor= element_blank(),
      strip.background= element_rect(fill = "gray95", color = "black", linewidth = 0.4),
      strip.text = element_text(face = "bold"),
      plot.margin  = margin(12, 12, 12, 12)
    )
}
theme_set(theme_publication())

save_pub <- function(p, filename, width = 12, height = 7, dpi = 600, make_pdf = TRUE) {
  png_file <- file.path(output_dir, paste0(filename, ".png"))
  pdf_file <- file.path(output_dir, paste0(filename, ".pdf"))
  
  ggsave(png_file, p, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", basename(png_file))
  
  if (isTRUE(make_pdf)) {
    ggsave(pdf_file, p, width = width, height = height, device = cairo_pdf, bg = "white")
    message("Saved: ", basename(pdf_file))
  }
}

# =============================================================================
# HELPERS: LOADING + STANDARDIZING
# =============================================================================
safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
}

find_metrics_file <- function(base_path, model_type = c("linear", "nn")) {
  model_type <- match.arg(model_type)
  
  candidates <- c(
    file.path(base_path, "aggregated_metrics.csv"),
    file.path(base_path, "all_genes_combined_metrics.csv"),
    file.path(base_path, "summary_statistics.csv"),
    file.path(base_path, "metrics_summary.csv"),
    file.path(base_path, "combined_metrics.csv"),
    file.path(base_path, "test_metrics.csv")
  )
  
  # DeepNN results sometimes nested
  if (model_type == "nn") {
    candidates <- c(
      candidates,
      Sys.glob(file.path(base_path, "*", "all_genes_combined_metrics.csv")),
      Sys.glob(file.path(base_path, "*", "*", "all_genes_combined_metrics.csv"))
    )
  }
  
  candidates <- unique(candidates)
  out <- candidates[file.exists(candidates)][1]
  if (is.na(out)) return(NULL)
  out
}

standardize_metrics_cols <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  # Common renames
  nm <- names(df)
  if (!"Model" %in% nm && "model" %in% nm) df <- df %>% rename(Model = model)
  if (!"Split" %in% nm && "split" %in% nm) df <- df %>% rename(Split = split)
  if (!"Gene" %in% nm && "gene" %in% nm) df <- df %>% rename(Gene = gene)
  
  # Harmonize common metric column variants
  if (!"R2" %in% names(df)) {
    if ("R_squared" %in% names(df)) df <- df %>% rename(R2 = R_squared)
    if ("Rsquared" %in% names(df)) df <- df %>% rename(R2 = Rsquared)
  }
  if (!"Spearman" %in% names(df)) {
    if ("spearman" %in% names(df)) df <- df %>% rename(Spearman = spearman)
    if ("rho" %in% names(df)) df <- df %>% rename(Spearman = rho)
  }
  if (!"RMSE" %in% names(df)) {
    if ("rmse" %in% names(df)) df <- df %>% rename(RMSE = rmse)
  }
  
  # Ensure Split exists
  if (!"Split" %in% names(df)) df$Split <- "Test"
  
  df
}

clean_model_names <- function(x) {
  x <- as.character(x)
  case_when(
    str_detect(x, "(?i)ols|ordinary") ~ "OLS",
    str_detect(x, "(?i)ridge")  ~ "Ridge",
    str_detect(x, "(?i)lasso")  ~ "Lasso",
    str_detect(x, "(?i)elastic")  ~ "ElasticNet",
    str_detect(x, "(?i)random|rf") ~ "RandomForest",
    str_detect(x, "(?i)deep|nn|neural") ~ "DeepNN",
    TRUE ~ x
  )
}

clip_metric_ranges <- function(df) {
  if ("Spearman" %in% names(df)) df <- df %>% mutate(Spearman = pmax(pmin(Spearman, 1), -1))
  if ("R2" %in% names(df))  df <- df %>% mutate(R2  = pmax(pmin(R2, 1), -1))
  df
}

load_metrics_folder <- function(dataset_name, base_path, model_type = c("linear", "nn")) {
  model_type <- match.arg(model_type)
  
  file <- find_metrics_file(base_path, model_type)
  if (is.null(file)) {
    message(" WARNING: no metrics file found in: ", base_path)
    return(NULL)
  }
  
  df <- safe_read_csv(file)
  df <- standardize_metrics_cols(df)
  
  if (is.null(df) || nrow(df) == 0) {
    message(" WARNING: could not read metrics: ", file)
    return(NULL)
  }
  
  df <- df %>%
    mutate(
      Dataset = dataset_name,
      Model = clean_model_names(Model)
    ) %>%
    filter(Split == "Test") %>%
    clip_metric_ranges()
  
  # Keep only expected models by folder type
  if (model_type == "linear") {
    df <- df %>% filter(Model %in% c("OLS", "Ridge", "Lasso", "ElasticNet", "RandomForest"))
  } else {
    df <- df %>% mutate(Model = "DeepNN") %>% filter(Model == "DeepNN")
  }
  
  message(" Loaded ", nrow(df), " rows from: ", file)
  df
}

load_dataset <- function(dataset_name, cfg) {
  message("\nLoading dataset: ", dataset_name)
  
  linear_df <- load_metrics_folder(dataset_name, cfg$linear_path, "linear")
  nn_df  <- load_metrics_folder(dataset_name, cfg$nn_path, "nn")
  
  out <- bind_rows(linear_df, nn_df)
  if (nrow(out) == 0) return(NULL)
  
  out %>%
    mutate(
      Dataset = factor(Dataset, levels = names(datasets)),
      Model  = factor(Model, levels = model_levels)
    )
}

# =============================================================================
# MAIN: LOAD ALL DATA
# =============================================================================
message("\n", strrep("=", 78))
message("STEP 085: Multi-Dataset Combined Performance Plots (Publication)")
message(strrep("=", 78), "\n")

all <- purrr::imap(datasets, ~ load_dataset(dataset_name = .y, cfg = .x)) %>% purrr::compact()

if (length(all) == 0) stop("ERROR: No datasets were successfully loaded.")

combined <- bind_rows(all)

message("\nCombined rows: ", nrow(combined))
message("Datasets: ", paste(levels(combined$Dataset), collapse = ", "))
message("Models:  ", paste(levels(combined$Model), collapse = ", "))

# Convenience summaries
median_by_group <- function(df, metric) {
  req <- c("Dataset", "Model", metric)
  if (!all(req %in% names(df))) return(tibble())
  df %>%
    filter(!is.na(.data[[metric]])) %>%
    group_by(Dataset, Model) %>%
    summarise(med = median(.data[[metric]], na.rm = TRUE), n = n(), .groups = "drop")
}

# =============================================================================
# FIG 1: SPEARMAN (Violin + slim box, no clutter)
# =============================================================================
if ("Spearman" %in% names(combined)) {
  sp <- combined %>% filter(!is.na(Spearman))
  if (nrow(sp) > 0) {
    
    sp_m <- median_by_group(sp, "Spearman")
    
    p1 <- ggplot(sp, aes(x = Model, y = Spearman, fill = Dataset)) +
      geom_violin(
        position = position_dodge(width = 0.82),
        alpha = 0.75, color = "black", linewidth = 0.25, trim = TRUE
      ) +
      geom_boxplot(
        position = position_dodge(width = 0.82),
        width = 0.18, outlier.shape = NA,
        alpha = 0.95, color = "black", linewidth = 0.35
      ) +
      geom_point(
        data = sp_m,
        aes(y = med, group = Dataset),
        position = position_dodge(width = 0.82),
        shape = 21, size = 2.1, stroke = 0.7,
        fill = "white", color = "black"
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
      scale_fill_manual(values = dataset_palette, name = "Dataset") +
      scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
      labs(
        title = "Test-set Spearman correlation by model and dataset",
        subtitle = "Distribution across genes (violin) with median (white point)",
        x = "Model",
        y = expression("Spearman correlation (" * rho * ")"),
        caption = "Boxplots show IQR; dashed line indicates 0 correlation."
      ) +
      theme(legend.position = "right") +
      theme(axis.text.x = element_text(face = "bold"),
            panel.grid.major.x = element_blank())
    
    save_pub(p1, "Figure1_Spearman_by_Model", width = 14, height = 8)
  }
}

# =============================================================================
# FIG 2: R2 (Box + light jitter)
# =============================================================================
if ("R2" %in% names(combined)) {
  r2 <- combined %>% filter(!is.na(R2))
  if (nrow(r2) > 0) {
    
    r2_m <- median_by_group(r2, "R2")
    
    p2 <- ggplot(r2, aes(x = Model, y = R2, fill = Dataset)) +
      geom_boxplot(
        position = position_dodge(width = 0.82),
        width = 0.55, outlier.shape = NA,
        alpha = 0.85, color = "black", linewidth = 0.35
      ) +
      geom_point(
        position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.82),
        alpha = 0.12, size = 0.7
      ) +
      geom_point(
        data = r2_m,
        aes(y = med, group = Dataset),
        position = position_dodge(width = 0.82),
        shape = 21, size = 2.1, stroke = 0.7,
        fill = "white", color = "black"
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
      scale_fill_manual(values = dataset_palette, name = "Dataset") +
      scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
      labs(
        title = expression("Test-set " * R^2 * " by model and dataset"),
        subtitle = "Boxplots with per-gene points (light jitter) and median (white point)",
        x = "Model",
        y = expression(R^2),
        caption = "Dashed line indicates 0 explained variance."
      ) +
      theme(legend.position = "right") +
      theme(axis.text.x = element_text(face = "bold"),
            panel.grid.major.x = element_blank())
    
    save_pub(p2, "Figure2_R2_by_Model", width = 14, height = 8)
  }
}

# =============================================================================
# FIG 3: RMSE (Box + robust trimming per group)
# =============================================================================
if ("RMSE" %in% names(combined)) {
  rmse <- combined %>% filter(!is.na(RMSE))
  if (nrow(rmse) > 0) {
    
    # Trim extreme tails within each (Dataset, Model): keep <= 99th percentile
    rmse_t <- rmse %>%
      group_by(Dataset, Model) %>%
      mutate(p99 = quantile(RMSE, 0.99, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(RMSE <= p99)
    
    rmse_m <- median_by_group(rmse_t, "RMSE")
    
    p3 <- ggplot(rmse_t, aes(x = Model, y = RMSE, fill = Dataset)) +
      geom_boxplot(
        position = position_dodge(width = 0.82),
        width = 0.55, outlier.shape = NA,
        alpha = 0.85, color = "black", linewidth = 0.35
      ) +
      geom_point(
        position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.82),
        alpha = 0.10, size = 0.7
      ) +
      geom_point(
        data = rmse_m,
        aes(y = med, group = Dataset),
        position = position_dodge(width = 0.82),
        shape = 21, size = 2.1, stroke = 0.7,
        fill = "white", color = "black"
      ) +
      scale_fill_manual(values = dataset_palette, name = "Dataset") +
      scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
      labs(
        title = "Test-set RMSE by model and dataset",
        subtitle = "Trimmed at the 99th percentile within each dataset-model group",
        x = "Model",
        y = "RMSE",
        caption = "Lower RMSE indicates better performance."
      ) +
      theme(legend.position = "right") +
      theme(axis.text.x = element_text(face = "bold"),
            panel.grid.major.x = element_blank())
    
    save_pub(p3, "Figure3_RMSE_by_Model", width = 14, height = 8)
  }
}

# =============================================================================
# FIG 4: FACET GRID (Dataset rows × Metric cols)
# =============================================================================
metrics_present <- intersect(c("R2", "RMSE", "Spearman"), names(combined))

if (length(metrics_present) > 0) {
  fd <- combined %>%
    pivot_longer(cols = all_of(metrics_present), names_to = "Metric", values_to = "Value") %>%
    filter(!is.na(Value)) %>%
    mutate(
      Metric = recode(Metric, R2 = "R²", RMSE = "RMSE", Spearman = "Spearman ρ"),
      Metric = factor(Metric, levels = c("R²", "Spearman ρ", "RMSE"))
    )
  
  # Use violin for Spearman, box for others
  p4 <- ggplot(fd, aes(x = Model, y = Value, fill = Model)) +
    geom_violin(
      data = fd %>% filter(Metric == "Spearman ρ"),
      alpha = 0.75, trim = FALSE, color = "black", linewidth = 0.25
    ) +
    geom_boxplot(
      data = fd %>% filter(Metric != "Spearman ρ"),
      alpha = 0.85, outlier.shape = NA, color = "black", linewidth = 0.35
    ) +
    facet_grid(Dataset ~ Metric, scales = "free_y") +
    labs(
      title = "Model performance across datasets and metrics",
      subtitle = "Rows = datasets; columns = metrics",
      x = "Model",
      y = "Metric value"
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
  
  save_pub(p4, "Figure4_Faceted_Performance", width = 16, height = 11)
}

# =============================================================================
# FIG 5: COMBINED R2 + SPEARMAN (two panels)
# =============================================================================
combo_metrics <- intersect(c("R2", "Spearman"), names(combined))
if (length(combo_metrics) > 0) {
  cd <- combined %>%
    pivot_longer(cols = all_of(combo_metrics), names_to = "Metric", values_to = "Value") %>%
    filter(!is.na(Value)) %>%
    mutate(
      Metric = recode(Metric, R2 = "R²", Spearman = "Spearman ρ"),
      Metric = factor(Metric, levels = c("R²", "Spearman ρ"))
    )
  
  p5 <- ggplot(cd, aes(x = Model, y = Value, fill = Dataset)) +
    geom_boxplot(
      position = position_dodge(width = 0.82),
      outlier.shape = NA, color = "black", linewidth = 0.35, alpha = 0.85
    ) +
    facet_wrap(~ Metric, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = dataset_palette, name = "Dataset") +
    labs(
      title = "Comparison of R² and Spearman correlation across datasets",
      x = "Model",
      y = "Metric value"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
  
  save_pub(p5, "Figure5_R2_and_Spearman", width = 14, height = 8)
}

# =============================================================================
# FIG 6: HEATMAP OF MEDIANS (direction-corrected so higher = better)
# =============================================================================
# For RMSE, invert direction by using -RMSE for display scale (still print true RMSE labels).
hm <- combined %>%
  group_by(Dataset, Model) %>%
  summarise(
    Median_R2 = if ("R2" %in% names(.)) median(R2, na.rm = TRUE) else NA_real_,
    Median_Spearman = if ("Spearman" %in% names(.)) median(Spearman, na.rm = TRUE) else NA_real_,
    Median_RMSE  = if ("RMSE" %in% names(.)) median(RMSE, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  pivot_longer(starts_with("Median_"), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = str_remove(Metric, "Median_"),
    Metric = recode(Metric, R2 = "R²", Spearman = "Spearman ρ", RMSE = "RMSE"),
    # direction-correct: bigger is better
    Value_for_fill = if_else(Metric == "RMSE", -Value, Value)
  )

if (nrow(hm) > 0) {
  p6 <- ggplot(hm, aes(x = Model, y = Dataset, fill = Value_for_fill)) +
    geom_tile(color = "white", linewidth = 0.9) +
    geom_text(aes(label = sprintf("%.2f", Value)), size = 4, fontface = "bold") +
    facet_wrap(~ Metric, ncol = 3, scales = "free") +
    scale_fill_gradient2(
      low = "#D55E00", mid = "white", high = "#0072B2",
      midpoint = 0, name = "Median\n(direction-corrected)"
    ) +
    labs(
      title = "Median performance summary (higher fill = better)",
      subtitle = "RMSE is direction-corrected (fill uses -RMSE), but labels show true RMSE",
      x = "Model",
      y = "Dataset"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray35"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  save_pub(p6, "Figure6_Median_Heatmap", width = 16, height = 7)
}

# =============================================================================
# SUMMARY TABLE (CSV)
# =============================================================================
summary_tbl <- combined %>%
  summarise(across(everything(), ~.x)) %>% invisible()

summary_out <- combined %>%
  pivot_longer(cols = intersect(c("R2", "Spearman", "RMSE"), names(combined)),
               names_to = "Metric", values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  group_by(Dataset, Model, Metric) %>%
  summarise(
    N_rows  = n(),
    N_genes = if ("Gene" %in% names(combined)) n_distinct(Gene) else NA_integer_,
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Dataset, Metric, desc(Median))

write_csv(summary_out, file.path(output_dir, "overall_summary_statistics.csv"))
message("\nSaved: overall_summary_statistics.csv")

# =============================================================================
# DONE
# =============================================================================
runtime <- difftime(Sys.time(), start_time, units = "mins")
message("\n", strrep("=", 78))
message("STEP 085 COMPLETE")
message("Output directory: ", output_dir)
message(sprintf("Runtime: %.2f minutes", as.numeric(runtime)))
message(strrep("=", 78), "\n")