library(tidyverse) 
library(data.table) # Crucial for reading large coefficient files fast 
library(ggpubr) 


# --- Custom Publication Theme ---
theme_pub <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(color = "grey30", size = base_size - 1),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "black", size = base_size - 1),
      legend.position = "top",
      legend.justification = "right",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 1)
    )
}

# --- Colors ---
cols <- c("RNA(Gene)" = "#377EB8", "ATAC(Peak)" = "#E41A1C")

# ============================================================================== 
# 1. CONFIGURATION: Define paths to your LINEAR/RF OUTPUT folders 
# ============================================================================== 
# Update BASE_RESULTS_DIR to point to your scMultiPreDICT output directory
# IMPORTANT: These folders must contain your 'coefficients.csv' or 'feature_importance.csv'

BASE_RESULTS_DIR <- "~/scMultiPreDICT_output/results/models"
DIMRED_METHOD <- "pca_lsi"       # "pca_lsi"
GENE_SET <- "HVG"                # Options: "HVG", "Random_genes"

# Define samples to analyze
SAMPLES <- c("E7.5_rep1", "E7.5_rep2", "T_Cells")
SAMPLE_LABELS <- c("E7.5_REP1", "E7.5_REP2", "T_Cells")

# Build dataset paths
datasets <- setNames(
  sapply(SAMPLES, function(s) file.path(BASE_RESULTS_DIR, "LINEAR_AND_TREE_BASED", s, DIMRED_METHOD, GENE_SET)),
  SAMPLE_LABELS
)

# Output directory for figures
OUTPUT_FIG_DIR <- file.path("~/scMultiPreDICT_output/paper_figures/feature_analysis", DIMRED_METHOD, GENE_SET)
dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)


load_dataset_coeffs <- function(dataset_name, base_path) {
  message(sprintf("Processing %s...", dataset_name))
  file_patterns <- c("coefficients_ols"="OLS", "coefficients_lasso"="Lasso", 
                     "coefficients_ridge"="Ridge", "coefficients_enet"="ElasticNet", 
                     "rf_importance"="RandomForest")
  
  all_data_list <- list()
  for (pattern in names(file_patterns)) {
    model_name <- file_patterns[[pattern]]
    files <- list.files(base_path, pattern = paste0(pattern, ".*\\.csv"), recursive = TRUE, full.names = TRUE)
    
    if (length(files) > 0) {
      model_df <- rbindlist(lapply(files, function(f) {
        dt <- fread(f, showProgress = FALSE)
        if ("Coefficient" %in% names(dt)) setnames(dt, "Coefficient", "Score")
        if ("Importance" %in% names(dt)) setnames(dt, "Importance", "Score")
        dt[, Target_Gene := basename(dirname(f))] 
        dt[, Model := model_name]
        return(dt[, .(Target_Gene, Feature, Score, Model)])
      }), fill = TRUE)
      all_data_list[[model_name]] <- model_df
    }
  }
  
  if (length(all_data_list) == 0) return(NULL)
  combined_dt <- rbindlist(all_data_list)
  combined_dt[, Dataset := dataset_name]
  combined_dt[, Feature_Type := ifelse(grepl(":", Feature) | grepl("^chr", Feature), "ATAC(Peak)", "RNA(Gene)")]
  return(combined_dt)
}

all_coeffs <- map_dfr(names(datasets), ~load_dataset_coeffs(., datasets[[.]]))

# ==============================================================================
# 3. CALCULATIONS
# ==============================================================================

# A. Calculate Denominator (Average Available Features)
dataset_averages <- all_coeffs %>%
  group_by(Dataset, Target_Gene, Feature_Type) %>%
  summarise(Total_Available = n_distinct(Feature), .groups = "drop") %>%
  group_by(Dataset, Feature_Type) %>%
  summarise(Avg_Input_Count = mean(Total_Available), .groups = "drop")

# B. Calculate Numerator (Selected in Top 30)
selected_counts <- all_coeffs %>%
  group_by(Dataset, Model, Target_Gene) %>%
  arrange(desc(abs(Score))) %>%
  slice_head(n = 30) %>% 
  ungroup() %>%
  count(Dataset, Model, Target_Gene, Feature_Type, name = "Selected_Count") %>%
  complete(nesting(Dataset, Model, Target_Gene), Feature_Type = c("RNA(Gene)", "ATAC(Peak)"), fill = list(Selected_Count = 0))

# C. Calculate Rates
plot_data <- selected_counts %>%
  left_join(dataset_averages, by = c("Dataset", "Feature_Type")) %>%
  mutate(
    # RENAMED METRIC: Selection Rate
    Selection_Rate = (Selected_Count / Avg_Input_Count) * 100
  )

# ==============================================================================
# 4. PLOT 1: RAW COUNTS (Separate PDF)
# ==============================================================================
stats_raw <- plot_data %>%
  group_by(Dataset, Model, Feature_Type) %>%
  summarise(Mean = mean(Selected_Count), SE = sd(Selected_Count)/sqrt(n()), .groups="drop")

p_raw <- ggplot(stats_raw, aes(x = Model, y = Mean, fill = Feature_Type)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(0.8), width = 0.25, linewidth = 0.5) +
  facet_wrap(~Dataset, nrow = 1) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Feature Composition (Raw Count)",
    subtitle = "Average number of features selected in the Top 30",
    y = "Count (Features)",
    x = NULL,
    fill = "Feature Type"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_FIG_DIR, "Fig_Feature_Raw_Counts.pdf"), p_raw, width = 12, height = 6)
message("Saved: Fig_Feature_Raw_Counts.pdf")

# ==============================================================================
# 5. PLOT 2: NORMALIZED SELECTION RATE (Separate PDF)
# ==============================================================================
stats_norm <- plot_data %>%
  group_by(Dataset, Model, Feature_Type) %>%
  summarise(Mean = mean(Selection_Rate), SE = sd(Selection_Rate)/sqrt(n()), .groups="drop")

p_norm <- ggplot(stats_norm, aes(x = Model, y = Mean, fill = Feature_Type)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(0.8), width = 0.25, linewidth = 0.5) +
  facet_wrap(~Dataset, nrow = 1) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Normalized Feature Selection",
    subtitle = "Rate of selection relative to available input features",
    y = "Selection Rate (% of Available)",
    x = NULL,
    fill = "Feature Type"
  ) +
  theme_pub() +
  theme(panel.border = element_blank(), axis.line = element_line(colour="black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_FIG_DIR, "Fig_Feature_Selection_Normalized.pdf"), p_norm, width = 12, height = 6)
message("Saved: Fig_Feature_Selection_Normalized.pdf")

print(p_norm)
# ==============================================================================
# 6. PLOT 3: INPUT DISTRIBUTION (Linear Scale, Separate PDF)
# ==============================================================================
dist_data <- all_coeffs %>%
  group_by(Dataset, Target_Gene, Feature_Type) %>%
  summarise(Feature_Count = n_distinct(Feature), .groups = "drop")

p_dist <- ggplot(dist_data, aes(x = Feature_Count, fill = Feature_Type, color = Feature_Type)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  facet_wrap(~Dataset, scales = "free") + # Allow x-axis to vary if needed
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  
  # LINEAR SCALE (No Log)
  scale_x_continuous(labels = scales::comma) +
  
  # Vertical Lines for Averages
  geom_vline(data = dataset_averages, aes(xintercept = Avg_Input_Count, color = Feature_Type),
             linetype = "dashed", linewidth = 0.8, alpha = 0.8) +
  
  labs(
    title = "Distribution of Input Features per Gene",
    subtitle = "Density of available features across target genes",
    y = "Density",
    x = "Number of Available Features",
    fill = "Feature Type"
  ) +
  theme_pub() +
  theme(
    legend.position = "top",
    # If the tail is very long, x-axis labels might overlap, so rotate them slightly
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
print(p_dist)
ggsave("Fig_Input_Feature_Distribution.pdf", p_dist, width = 14, height = 6)
message("Saved: Fig_Input_Feature_Distribution.pdf")