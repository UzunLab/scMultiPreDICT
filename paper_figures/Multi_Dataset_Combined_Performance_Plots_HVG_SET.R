#!/usr/bin/env Rscript
# ============================================================================
# Multi-Dataset Combined Performance Plots (Publication Ready)
# ============================================================================

start_time <- Sys.time()

# ============================================================================
# 1. LOAD LIBRARIES
# ============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(viridis)
  library(ggrepel)
  library(scales)
  library(ggpubr)    # For statistical significance
  library(gghalves)  # For Raincloud plots (Half-violins)
  library(cowplot)   # For arranging complex grids
  library(Cairo)
  library(gridExtra)
})

# ============================================================================
# 2. CONFIGURATION - EDIT THIS SECTION
# ============================================================================
# Update BASE_RESULTS_DIR to point to your scMultiPreDICT output directory

BASE_RESULTS_DIR <- "~/scMultiPreDICT_output/results/models"
DIMRED_METHOD <- "pca_lsi"  # Options: "pca_lsi", "wnn", "scvi_peakvi", "multivi"
GENE_SET <- "HVG"           # Options: "HVG", "Random_genes"

# Define samples to analyze
SAMPLES <- c("E7.5_rep1", "E7.5_rep2", "T_Cells")
SAMPLE_LABELS <- c("E7.5_REP1", "E7.5_REP2", "T_Cells")

# Build dataset paths
datasets <- setNames(lapply(SAMPLES, function(s) {
  list(
    linear_path = file.path(BASE_RESULTS_DIR, "LINEAR_AND_TREE_BASED", s, DIMRED_METHOD, GENE_SET),
    nn_path = file.path(BASE_RESULTS_DIR, "NEURAL_NETWORKS", s, DIMRED_METHOD, "Three_hidden_layer", "Three_hidden_layer", GENE_SET)
  )
}), SAMPLE_LABELS)

# Output directory
output_dir <- file.path("~/scMultiPreDICT_output/paper_figures", DIMRED_METHOD, GENE_SET)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 3. FIGURES THEME & PALETTES
# ============================================================================

# Dataset palette 
dataset_cols <- c("E7.5_REP1" = "#E69F00", "E7.5_REP2" = "#56B4E9", "T_Cells" = "#009E73")
# Ash color used for jitter points ('glittering')
ash_col <- "#3F3837"
model_cols   <- c("OLS" = "#999999", "Ridge" = "#0072B2", "Lasso" = "#009E73", 
                  "ElasticNet" = "#E69F00", "RandomForest" = "#D55E00", "DeepNN" = "#CC79A7")

theme_nature <- function(base_size = 12, base_family = "Arial") {
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
theme_set(theme_nature())

# ============================================================================
# 4. DATA LOADING FUNCTIONS
# ============================================================================

find_file <- function(base_path) {
  # Priorities: specific combined metrics -> generic summary -> recursive search
  patterns <- c("all_genes_combined_metrics.csv", "metrics_summary.csv", "combined_metrics.csv", "aggregated_metrics.csv")
  for (p in patterns) {
    f <- file.path(base_path, p)
    if (file.exists(f)) return(f)
  }
  # Recursive search if standard path fails
  fs <- list.files(base_path, pattern = "combined_metrics.csv", recursive = TRUE, full.names = TRUE)
  if (length(fs) > 0) return(fs[1])
  return(NULL)
}

load_data_safe <- function(dataset_name, config) {
  # Load Linear
  f_lin <- find_file(config$linear_path)
  d_lin <- if (!is.null(f_lin)) read_csv(f_lin, show_col_types = FALSE) else NULL
  
  # Load NN
  f_nn <- find_file(config$nn_path)
  d_nn <- if (!is.null(f_nn)) read_csv(f_nn, show_col_types = FALSE) else NULL
  
  # Combine
  if (is.null(d_lin) & is.null(d_nn)) return(NULL)
  combined <- bind_rows(d_lin, d_nn) %>%
    mutate(Dataset = dataset_name) %>%
    filter(Split == "Test") # Ensure we only use Test data
  
  # Clean Model Names
  combined <- combined %>%
    rename_with(str_to_title) %>% # Fix column case issues (model -> Model)
    # Ensure consistent metric column names (e.g. RMSE -> RMSE not Rmse)
    rename_with(~ ifelse(. == "Rmse", "RMSE", .)) %>%
    mutate(Model = case_when(
      str_detect(Model, "(?i)ols") ~ "OLS",
      str_detect(Model, "(?i)ridge") ~ "Ridge",
      str_detect(Model, "(?i)lasso") ~ "Lasso",
      str_detect(Model, "(?i)elastic") ~ "ElasticNet",
      str_detect(Model, "(?i)random|rf") ~ "RandomForest",
      str_detect(Model, "(?i)deep|nn") ~ "DeepNN",
      TRUE ~ Model
    )) %>%
    filter(Model %in% names(model_cols)) # Filter only known models
  
  return(combined)
}

# ============================================================================
# 5. MAIN EXECUTION
# ============================================================================

cat("Loading Datasets...\n")
all_data <- map_dfr(names(datasets), ~load_data_safe(., datasets[[.]]))

# Factor Ordering (Crucial for plotting order)
all_data$Model <- factor(all_data$Model, levels = c("OLS", "Lasso", "Ridge", "ElasticNet", "RandomForest", "DeepNN"))
all_data$Dataset <- factor(all_data$Dataset, levels = names(datasets))

cat(sprintf("Loaded %d rows across %d datasets.\n", nrow(all_data), n_distinct(all_data$Dataset)))

# ============================================================================
# 6. GENERATE FIGURES
# ============================================================================

# --- FIGURE A: Raincloud Plot for Spearman Correlation ---
# This is better than a simple violin because it shows raw data points + density + boxplot
plot_data_rain <- all_data
if (nrow(plot_data_rain) > 5000) {
  set.seed(42)
  plot_data_rain <- plot_data_rain %>% group_by(Dataset, Model) %>% sample_n(size = pmin(500, n()), replace = FALSE) %>% ungroup()
}

p_raincloud <- ggplot(plot_data_rain, aes(x = Model, y = Spearman, fill = Dataset)) +
  geom_half_violin(
    position = position_nudge(x = 0.2, y = 0),
    side = "r", adjust = 1.2, trim = FALSE, alpha = 0.7, color = NA
  ) +
  geom_boxplot(
    width = 0.35, outlier.shape = NA, alpha = 0.95,
    position = position_nudge(x = 0.2), color = "black", size = 0.25
  ) +
  # black median line for clarity (shorter than full box)
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 0, size = 0.6, color = "black", position = position_nudge(x = 0.2)) +
  geom_jitter(shape = 21, fill = ash_col, color = "black", stroke = 0.08, position = position_jitter(width = 0.06, height = 0), size = 0.6, alpha = 0.25, show.legend = FALSE) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.2, color = "black", fill = "white", position = position_nudge(x = 0.2)) +
  facet_wrap(~Dataset, nrow = 1) +
  scale_fill_manual(values = dataset_cols, name = "Dataset") +
  scale_color_manual(values = dataset_cols) +
  scale_y_continuous(limits = c(-0.1, 1.0), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Spearman Correlation Distribution (HVGs)",
    y = "Spearman Correlation",
    x = NULL
  ) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "Fig4_WNN_Correlation_Combined_Modality_hvg.pdf"), p_raincloud, width = 12, height = 7, device = cairo_pdf)
ggsave(file.path(output_dir, "Fig4_WNN_Correlation_Combined_Modality_hvg.png"), p_raincloud, width = 12, height = 7, dpi = 600)


# --- FIGURE B: Boxplots for R² and RMSE 

# Ensure metrics exist and compute simple summary table for heatmaps
metrics_present <- c("R2", "Spearman")
if ("RMSE" %in% colnames(all_data)) metrics_present <- c(metrics_present, "RMSE")

summary_stats <- all_data %>%
  group_by(Dataset, Model) %>%
  summarise(across(all_of(metrics_present), ~ mean(.x, na.rm = TRUE), .names = "Mean_{.col}"), .groups = "drop")

# Boxplot: R2
if ("R2" %in% colnames(all_data)) {
  p_box_r2 <- ggplot(all_data, aes(x = Model, y = R2, fill = Dataset)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6) +
      stat_summary(fun = median, geom = "point", shape = 21, size = 2.2, fill = "white", color = "black") +
      geom_jitter(shape = 21, fill = ash_col, color = "black", stroke = 0.08, width = 0.12, size = 0.6, alpha = 0.25, show.legend = FALSE) +
    facet_wrap(~Dataset, nrow = 1, scales = "fixed") +
    scale_fill_manual(values = dataset_cols, name = "Dataset") +
    scale_color_manual(values = dataset_cols) +
    # guides(fill = guide_legend(title = "Dataset")) +
    labs(title = "R-Squared Distribution (HVGs)", y = "R²", x = NULL) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(-0.5, 1.0), breaks = seq(-0.5, 1.0, 0.25))

  # mean labels removed for cleaner appearance

  ggsave(file.path(output_dir, "Fig4_WNN_R2_Combined_Modality_hvg.pdf"), p_box_r2, width = 12, height = 8, device = cairo_pdf)
  ggsave(file.path(output_dir, "Fig4_WNN_R2_Combined_Modality_hvg.png"), p_box_r2, width = 12, height = 8, dpi = 600)
}

# Boxplot: RMSE
if ("RMSE" %in% colnames(all_data)) {
  p_box_rmse <- ggplot(all_data, aes(x = Model, y = RMSE, fill = Dataset)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6) +
      stat_summary(fun = median, geom = "point", shape = 21, size = 2.2, fill = "white", color = "black") +
      geom_jitter(shape = 21, fill = ash_col, color = "black", stroke = 0.08, width = 0.12, size = 0.6, alpha = 0.25, show.legend = FALSE) +
    facet_wrap(~Dataset, nrow = 1, scales = "fixed") +
    scale_fill_manual(values = dataset_cols, name = "Dataset") +
    scale_color_manual(values = dataset_cols) +
    labs(title = "Root Mean Squared Error (HVGs)", y = "RMSE", x = NULL) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

  # mean labels removed for cleaner appearance

  ggsave(file.path(output_dir, "Fig4_WNN_RMSE_Combined_Modality_hvg.pdf"), p_box_rmse, width = 12, height = 8, device = cairo_pdf)
  ggsave(file.path(output_dir, "Fig4_WNN_RMSE_Combined_Modality_hvg.png"), p_box_rmse, width = 12, height = 8, dpi = 600)
}

# --- FIGURE C: Comparison Heatmaps for all metrics ---

# Prepare long-form mean values per Metric/Dataset/Model
heatmap_dat <- all_data %>%
  pivot_longer(cols = all_of(metrics_present), names_to = "Metric", values_to = "Value") %>%
  group_by(Metric, Dataset, Model) %>%
  summarise(Median = median(Value, na.rm = TRUE), .groups = "drop") %>%
  group_by(Metric) %>%
  mutate(Thresh = (max(Median, na.rm = TRUE) + min(Median, na.rm = TRUE)) / 2,
         TextColor = ifelse(is.na(Median), "black", ifelse(Median >= Thresh, "white", "black"))) %>%
  ungroup()

p_heat_all <- ggplot(heatmap_dat, aes(x = Dataset, y = Model, fill = Median)) +
  geom_tile(color = "white", height = 0.92) +
  geom_text(aes(label = sprintf("%.2f", Median), color = TextColor), fontface = "bold", size = 4.0, show.legend = FALSE) +
  facet_wrap(~Metric, nrow = 1, scales = "free_x") +
  scale_fill_viridis_c(option = "maco", direction = 1, na.value = "grey80", name = "Median") +
  scale_color_identity() +
  labs(title = "Model Performance Overview (Median Metrics)", x = NULL, y = NULL) +
  theme_nature(base_size = 12) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "Fig3_Heatmap_AllMetrics.pdf"), p_heat_all, device = cairo_pdf,
       width = 10, height = 5 + length(metrics_present) * 0.3, units = "in")
ggsave(file.path(output_dir, "Fig3_Heatmap_AllMetrics.png"), p_heat_all, width = 10, height = 5 + length(metrics_present) * 0.3, units = "in", dpi = 600)

# Also save a separate heatmap file per metric for easier figure assembly
for (m in unique(heatmap_dat$Metric)) {
  sub_dat <- heatmap_dat %>% filter(Metric == m)
  p_m <- ggplot(sub_dat, aes(x = Dataset, y = Model, fill = Median)) +
    geom_tile(color = "white", height = 0.92) +
    geom_text(aes(label = sprintf("%.2f", Median), color = TextColor), fontface = "bold", size = 4.0, show.legend = FALSE) +
    scale_fill_viridis_c(option = "maco", direction = 1, na.value = "grey80", name = "Median") +
    scale_color_identity() +
    labs(title = paste0("Median ", m), x = NULL, y = NULL) +
    theme_nature(base_size = 12) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  
  fname_base <- paste0("Fig3_Heatmap_", gsub("[^A-Za-z0-9]", "", m))
  ggsave(file.path(output_dir, paste0(fname_base, ".pdf")), p_m, device = cairo_pdf,
         width = 5, height = 4, units = "in")
  ggsave(file.path(output_dir, paste0(fname_base, ".png")), p_m, width = 5, height = 4, units = "in", dpi = 600)
}


# --- FIGURE C: Comparison Heatmap ---
# 

heatmap_dat2 <- summary_stats %>%
  select(Dataset, Model, Median_R2) %>%
  mutate(Model = factor(Model, levels = levels(all_data$Model)), Dataset = factor(Dataset, levels = levels(all_data$Dataset)))

p_heat <- ggplot(heatmap_dat2, aes(x = Dataset, y = Model, fill = Median_R2)) +
  geom_tile(color = "white", height = 0.92) +
  geom_text(aes(label = sprintf("%.2f", Median_R2)), color = "white", fontface = "bold", size = 4) +
  scale_x_discrete(position = "bottom", labels = dataset_labels) +
  scale_fill_viridis_c(option = "maco", name = "Median R²", direction = 1, na.value = "grey80") +
  labs(title = "Model Performance Overview (Median R²)", x = NULL, y = NULL) +
  theme_nature(base_size = 11) +
  theme(
    axis.text = element_text(size = 9, color = "black", angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )
ggsave(file.path(output_dir, "Fig3_Heatmap_Summary_Correlation.pdf"), p_heat, device = cairo_pdf,
       width = 5, height = 3, units = "in", family = "Arial")
embed_pdf(file.path(output_dir, "Fig3_Heatmap_Summary_Correlation.pdf"))
ggsave(file.path(output_dir, "Fig3_Heatmap_Summary_Correlation.png"), p_heat, width = 6, height = 5, units = "in", dpi = 600)

# --- COMPLETE ---
cat("Figures saved to:", output_dir, "\n")
