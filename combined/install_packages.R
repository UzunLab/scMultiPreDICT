#!/usr/bin/env Rscript
# ============================================================================
# scMultiPreDICT - R Package Installation Script
# ============================================================================
#
# Description:
#   Automated installation of all R package dependencies required for the
#   scMultiPreDICT analysis pipeline.
#
# Usage:
#   source("install_packages.R")  # From R console
#   Rscript install_packages.R    # From command line
#
# ============================================================================

cat("============================================================\n")
cat("      scMultiPreDICT R Package Installation                 \n")
cat("============================================================\n\n")

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
  "RColorBrewer",
  "cowplot",
  "reshape2",
  
  # Single-cell analysis
  
  # Machine learning
  "glmnet",
  "ranger",
  
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
# Seurat (CRAN first, GitHub fallback)
# ============================================================================

install_matrix <- function() {
  if (requireNamespace("Matrix", quietly = TRUE)) {
    cat("  Already installed: Matrix\n")
    return(invisible(TRUE))
  }

  cat("  Installing: Matrix (CRAN)\n")
  tryCatch({
    install.packages("Matrix", repos = "https://cloud.r-project.org", quiet = TRUE)
  }, error = function(e) {
    cat(sprintf("  CRAN install error for Matrix: %s\n", e$message))
  })

  if (requireNamespace("Matrix", quietly = TRUE)) {
    cat("  Matrix installed successfully.\n")
  } else {
    cat("  Matrix installation failed.\n")
  }
}

install_mass <- function() {
  if (requireNamespace("MASS", quietly = TRUE)) {
    cat("  Already installed: MASS\n")
    return(invisible(TRUE))
  }

  cat("  Installing: MASS (CRAN)\n")
  tryCatch({
    install.packages("MASS", repos = "https://cloud.r-project.org", quiet = TRUE)
  }, error = function(e) {
    cat(sprintf("  CRAN install error for MASS: %s\n", e$message))
  })

  if (requireNamespace("MASS", quietly = TRUE)) {
    cat("  MASS installed successfully.\n")
  } else {
    cat("  MASS installation failed.\n")
  }
}

install_caret <- function() {
  if (requireNamespace("caret", quietly = TRUE)) {
    cat("  Already installed: caret\n")
    return(invisible(TRUE))
  }

  cat("  Installing: caret (CRAN)\n")
  tryCatch({
    install.packages("caret", repos = "https://cloud.r-project.org", quiet = TRUE)
  }, error = function(e) {
    cat(sprintf("  CRAN install error for caret: %s\n", e$message))
  })

  if (requireNamespace("caret", quietly = TRUE)) {
    cat("  caret installed successfully.\n")
  } else {
    cat("  caret installation failed.\n")
  }
}

install_seurat <- function() {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    cat("  Already installed: Seurat\n")
    return(invisible(TRUE))
  }

  cat("  Installing: Seurat (CRAN)\n")
  tryCatch({
    install.packages("Seurat", repos = "https://cloud.r-project.org", quiet = TRUE)
  }, error = function(e) {
    cat(sprintf("  CRAN install error for Seurat: %s\n", e$message))
  })

  if (requireNamespace("Seurat", quietly = TRUE)) {
    cat("  Seurat installed successfully from CRAN.\n")
    return(invisible(TRUE))
  }

  cat("  CRAN install failed or incomplete. Trying GitHub branch 'seurat5'...\n")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }

  tryCatch({
    remotes::install_github("satijalab/seurat", ref = "seurat5", quiet = TRUE, upgrade = "never")
  }, error = function(e) {
    cat(sprintf("  GitHub install error for Seurat: %s\n", e$message))
  })

  if (requireNamespace("Seurat", quietly = TRUE)) {
    cat("  Seurat installed successfully from GitHub.\n")
  } else {
    cat("  Seurat installation failed from both CRAN and GitHub.\n")
  }
}

cat("\nInstalling Seurat...\n")
cat("Ensuring Matrix is installed first...\n")
install_matrix()
cat("Ensuring MASS is installed for caret dependencies...\n")
install_mass()
cat("Installing caret after Matrix/MASS...\n")
install_caret()
install_seurat()

# ============================================================================
# Bioconductor Packages
# ============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
}

bioc_packages <- c(
  "Signac",
  "GenomicRanges",
  "GenomeInfoDb",
  "IRanges",
  "S4Vectors"
)

install_rtracklayer <- function() {
  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    cat("  Already installed: rtracklayer\n")
    return(invisible(TRUE))
  }

  cat("  Installing: rtracklayer (conda: bioconda::bioconductor-rtracklayer)\n")
  conda_bin <- Sys.which("conda")
  conda_ok <- FALSE

  if (nzchar(conda_bin)) {
    conda_status <- tryCatch({
      system2(conda_bin, args = c("install", "-y", "bioconda::bioconductor-rtracklayer"))
    }, error = function(e) {
      cat(sprintf("  Conda install error for rtracklayer: %s\n", e$message))
      1L
    })
    conda_ok <- identical(conda_status, 0L)
  } else {
    cat("  Conda not found in PATH; skipping conda install path.\n")
  }

  if (conda_ok && requireNamespace("rtracklayer", quietly = TRUE)) {
    cat("  rtracklayer installed successfully via conda.\n")
    return(invisible(TRUE))
  }

  cat("  Conda install failed or package still unavailable. Trying BiocManager...\n")
  tryCatch({
    BiocManager::install("rtracklayer", update = FALSE, ask = FALSE)
  }, error = function(e) {
    cat(sprintf("  BiocManager install error for rtracklayer: %s\n", e$message))
  })

  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    cat("  rtracklayer installed successfully via BiocManager.\n")
  } else {
    cat("  rtracklayer installation failed from both conda and BiocManager.\n")
  }
}

cat("\nInstalling Bioconductor packages...\n")
install_rtracklayer()
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing: %s\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("  Already installed: %s\n", pkg))
  }
}

# ============================================================================
# Species-Specific Annotation Packages
# ============================================================================

cat("\n============================================================\n")
cat("Species-Specific Annotation Packages\n")
cat("============================================================\n")
cat("Install ONE of the following based on your organism:\n\n")

cat("For Mus musculus (mouse):\n")
cat("  BiocManager::install('EnsDb.Mmusculus.v79')\n")
cat("  BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')\n\n")

cat("For Homo sapiens (human):\n")
cat("  BiocManager::install('EnsDb.Hsapiens.v86')\n")
cat("  BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')\n\n")

# Interactive species selection
species_choice <- readline(prompt = "Install annotation packages? (mouse/human/skip): ")

if (tolower(species_choice) == "mouse") {
  cat("Installing mouse annotation packages...\n")
  BiocManager::install("EnsDb.Mmusculus.v79", update = FALSE, ask = FALSE)
  BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", update = FALSE, ask = FALSE)
} else if (tolower(species_choice) == "human") {
  cat("Installing human annotation packages...\n")
  BiocManager::install("EnsDb.Hsapiens.v86", update = FALSE, ask = FALSE)
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update = FALSE, ask = FALSE)
} else {
  cat("Skipping species-specific packages. Install manually when needed.\n")
}

# ============================================================================
# Installation Verification
# ============================================================================

cat("\n============================================================\n")
cat("Installation Verification\n")
cat("============================================================\n")

all_packages <- c("Seurat", "Matrix", "MASS", "caret", "rtracklayer", cran_packages, bioc_packages)
missing <- c()

for (pkg in all_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing <- c(missing, pkg)
  }
}

if (length(missing) > 0) {
  cat("\nWARNING: The following packages failed to install:\n")
  cat(paste("  -", missing, collapse = "\n"), "\n")
  cat("\nPlease install these packages manually.\n")
} else {
  cat("\nAll core packages installed successfully.\n")
}

# ============================================================================
# Session Information
# ============================================================================

cat("\n============================================================\n")
cat("Session Information\n")
cat("============================================================\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform: %s\n", R.version$platform))
if (requireNamespace("BiocManager", quietly = TRUE)) {
  cat(sprintf("BiocManager version: %s\n", as.character(packageVersion("BiocManager"))))
}
if (requireNamespace("Seurat", quietly = TRUE)) {
  cat(sprintf("Seurat version: %s\n", as.character(packageVersion("Seurat"))))
}
if (requireNamespace("Signac", quietly = TRUE)) {
  cat(sprintf("Signac version: %s\n", as.character(packageVersion("Signac"))))
}

cat("\n============================================================\n")
cat("Installation Complete\n")
cat("============================================================\n")
