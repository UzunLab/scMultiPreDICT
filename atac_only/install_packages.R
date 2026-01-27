# scMultiPreDICT R Package Requirements
# ============================================================================
# Run this script to install all required R packages
#
# Usage:
#   Interactive: source("install_packages.R")
#   Command line: Rscript install_packages.R --species=mouse
#                 Rscript install_packages.R --species=human
#                 Rscript install_packages.R  (skips species packages)
# ============================================================================

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
species_arg <- NULL
for (arg in args) {
  if (grepl("^--species=", arg)) {
    species_arg <- sub("^--species=", "", arg)
  }
}

cat("============================================\n")
cat("Installing scMultiPreDICT R Dependencies\n")
cat("============================================\n\n")

# ============================================================================
# CRAN Packages
# ============================================================================
cran_packages <- c(
  # Data manipulation
  "dplyr",
  "tidyr",
  "readr",
  "purrr",
  "stringr",
  "tibble",
  "forcats",
  
  # Visualization
  "ggplot2",
  "patchwork",
  "viridis",
  "scales",
  "ggrepel",
  "Cairo",
  "RColorBrewer",
  "cowplot",
  "reshape2",
  
  # Single-cell analysis
  "Seurat",
  "Signac",
  "Matrix",
  
  # Machine learning
  "glmnet",
  "ranger",
  "caret",
  
  # Utilities
  "RANN",
  "irlba",
  "reticulate",
  "parallel",
  "doParallel",
  "foreach",
  "matrixStats"
)

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(sprintf("  Already installed: %s\n", pkg))
  }
}

# ============================================================================
# Bioconductor Packages
# ============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "GenomicRanges",
  "GenomeInfoDb",
  "rtracklayer",
  "IRanges",
  "S4Vectors"
)

cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing: %s\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("  Already installed: %s\n", pkg))
  }
}

# ============================================================================
# Species-Specific Packages
# ============================================================================
cat("\n============================================\n")
cat("Species-Specific Packages\n")
cat("============================================\n")
cat("Install ONE of the following based on your data:\n\n")

cat("For MOUSE data:\n")
cat("  BiocManager::install('EnsDb.Mmusculus.v79')\n")
cat("  BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')\n\n")

cat("For HUMAN data:\n")
cat("  BiocManager::install('EnsDb.Hsapiens.v86')\n")
cat("  BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')\n\n")

# Determine species choice (command-line arg or interactive prompt)
if (!is.null(species_arg)) {
  species_choice <- species_arg
  cat(sprintf("Species specified via command line: %s\n", species_choice))
} else if (interactive()) {
  species_choice <- readline(prompt = "Install species packages? (mouse/human/skip): ")
} else {
  cat("No --species argument provided. Skipping species packages.\n")
  cat("To install, run: Rscript install_packages.R --species=mouse\n")
  cat("             or: Rscript install_packages.R --species=human\n")
  species_choice <- "skip"
}

if (tolower(species_choice) == "mouse") {
  cat("Installing mouse annotation packages...\n")
  tryCatch({
    BiocManager::install("EnsDb.Mmusculus.v79", update = FALSE, ask = FALSE)
    BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", update = FALSE, ask = FALSE)
  }, error = function(e) {
    cat(sprintf("WARNING: Failed to install mouse packages: %s\n", e$message))
  })
} else if (tolower(species_choice) == "human") {
  cat("Installing human annotation packages...\n")
  tryCatch({
    BiocManager::install("EnsDb.Hsapiens.v86", update = FALSE, ask = FALSE)
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update = FALSE, ask = FALSE)
  }, error = function(e) {
    cat(sprintf("WARNING: Failed to install human packages: %s\n", e$message))
  })
} else {
  cat("Skipping species packages. Install manually when needed.\n")
}

# ============================================================================
# Verification
# ============================================================================
cat("\n============================================\n")
cat("Verifying Installation\n")
cat("============================================\n")

all_packages <- c(cran_packages, bioc_packages)
missing <- c()

for (pkg in all_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing <- c(missing, pkg)
  }
}

if (length(missing) > 0) {
  cat("\nWARNING: The following packages failed to install:\n")
  cat(paste("  -", missing, collapse = "\n"), "\n")
  cat("\nTry installing them manually.\n")
} else {
  cat("\n✓ All core packages installed successfully!\n")
}

# ============================================================================
# Session Info
# ============================================================================
cat("\n============================================\n")
cat("Session Information\n")
cat("============================================\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform: %s\n", R.version$platform))

# Safely print package versions
for (pkg in c("BiocManager", "Seurat", "Signac")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("%s version: %s\n", pkg, as.character(packageVersion(pkg))))
  } else {
    cat(sprintf("%s: NOT INSTALLED\n", pkg))
  }
}

cat("\n============================================\n")
if (length(missing) > 0) {
  cat("Installation Complete (with warnings)\n")
} else {
  cat("Installation Complete!\n")
}
cat("============================================\n")
