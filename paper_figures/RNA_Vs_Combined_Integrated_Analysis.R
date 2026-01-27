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
# This script compares RNA-only vs Combined (RNA+ATAC) performance
# Update BASE_RESULTS_DIR to point to your scMultiPreDICT output directory

RNA_RESULTS_DIR <- "~/scMultiPreDICT_output/rna_only/results/models"
RNA_METHOD <- "pca"              
GENE_SET <- "HVG"                # Options: "HVG", "Random_genes"

# Define samples to analyze
SAMPLES <- c("E7.5_rep1", "E7.5_rep2", "T_Cells")
SAMPLE_LABELS <- c("E7.5_REP1", "E7.5_REP2", "T_Cells")

# Build RNA-only dataset paths
pca_rna <- setNames(lapply(SAMPLES, function(s) {
  list(
    linear_path = file.path(RNA_RESULTS_DIR, "LINEAR_AND_TREE_BASED", s, RNA_METHOD, GENE_SET),
    nn_path = file.path(RNA_RESULTS_DIR, "NEURAL_NETWORKS", s, RNA_METHOD, "Three_hidden_layer", GENE_SET)
  )
}), SAMPLE_LABELS)

# Output directory
output_dir <- "~/scMultiPreDICT_output/paper_figures/rna_vs_combined/"
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

# PCA + LSI
cat("Loading Datasets...\n")
pca_rna_data <- map_dfr(names(pca_rna), ~load_data_safe(., pca_rna[[.]]))

# Factor Ordering (Crucial for plotting order)
pca_rna_data$Model <- factor(pca_rna_data$Model, levels = c("OLS", "Lasso", "Ridge", "ElasticNet", "RandomForest", "DeepNN"))
pca_rna_data$Dataset <- factor(pca_rna_data$Dataset, levels = names(pca_rna))

cat(sprintf("Loaded %d rows across %d datasets.\n", nrow(pca_rna_data), n_distinct(pca_rna_data$Dataset)))


# ============================================================================
# 6. Integrate all the dataset into a big dataframe
# ============================================================================
# Add an Integration column to each dimensionality strategy

# Select columns needed for plotting
df_main_rna <- pca_rna_data %>% 
  filter(
    Split == "Test",
    Geneset == "HVG"
  )

# Collapse across models
df_gene_level_rna <- df_main_rna %>% 
  group_by(Dataset, Gene) %>% 
  summarise(
    rna_spearman = median(Spearman, na.rm = TRUE),
    .groups = "drop"
  )

# Pool across datasets
df_delta <- 
  df_gene_level %>% 
  left_join(df_gene_level_rna, by = c("Dataset", "Gene")) %>% 
  mutate(delta_spear = spearman - rna_spearman)


# ============================================================================
# 6. GENERATE FIGURES
# ============================================================================

pC <- ggplot(df_delta, aes(x=Integrated_Analysis, y=delta_spear, fill = Integrated_Analysis)) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  geom_violin(trim=FALSE, alpha=0.75) + 
  geom_boxplot(width=0.18, outlier.shape=NA) + 
  labs(title = "Multimodal vs RNA-only Prediction", x=NULL, y="Δ Spearman (Multimodal − RNA-only)") + 
  theme(axis.text.x = element_text(angle=35, hjust=1)) 
pC
ggsave(file.path(output_dir, "Fig4_Multimodal_vs_RNA_only.pdf"), pC, device = cairo_pdf,
       width = 10, height = 5 + length(metrics_present) * 0.3, units = "in")
ggsave(file.path(output_dir, "Fig4_Multimodal_vs_RNA_only.png"), pC, width = 10, height = 5 + length(metrics_present) * 0.3, units = "in", dpi = 600)


# Percentage of Genes Imporved if any
df_improve <- df_delta %>% 
  mutate(status = case_when( 
    delta_spear >  0.01 ~ "Improved", 
    delta_spear < -0.01 ~ "Worse", 
    TRUE ~ "No change" 
  )) %>% 
  count(Integrated_Analysis, status) %>% 
  group_by(Integrated_Analysis) %>% 
  mutate(percentage = (n / sum(n)) * 100,
         percentage_label = round(percentage, 1)) 

pE <- ggplot(df_improve, aes(x=Integrated_Analysis, y=percentage, fill=status)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black") + 
  geom_text(aes(label = paste0(percentage_label, "%")),
            position = position_stack(vjust = 0.5), size = 4, color = "white") +
  scale_fill_manual(values = c("Improved" = "#E41A1C",
                               "Worse" = "#377EB8",
                               "No change" = "#999999")) +
  labs(title = "Impact of Combined Modalities Per Target Gene", x=NULL, y="Percentage of Genes (%)") + 
  theme_nature() 

pE 

ggsave(file.path(output_dir, "Percentage_of_Genes_Improved.pdf"), pE, device = cairo_pdf,
       width = 7, height = 4, units = "in")
ggsave(file.path(output_dir, "Percentage_of_Genes_Improved.png"), pE, width = 7, height = 4, units = "in", dpi = 600)
