#!/usr/bin/env Rscript
# ============================================================================
# scMultiPreDICT - Installation Test Script
# ============================================================================
#
# Description:
#   This script verifies that all required R packages are installed and
#   checks the Python environment for optional autoencoder-based methods.
#
# Usage:
#   Rscript test_installation.R
#
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("         scMultiPreDICT Installation Test                   \n")
cat("============================================================\n\n")

# Track test results
tests_passed <- 0
tests_failed <- 0
warnings_list <- c()

# Helper function to test package availability
test_package <- function(pkg, required = TRUE) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- as.character(packageVersion(pkg))
    cat(sprintf("  [✓] %s (%s)\n", pkg, version))
    return(TRUE)
  } else {
    if (required) {
      cat(sprintf("  [✗] %s - NOT FOUND (REQUIRED)\n", pkg))
    } else {
      cat(sprintf("  [!] %s - NOT FOUND (optional)\n", pkg))
    }
    return(FALSE)
  }
}

# ============================================================================
# Test Core R Packages
# ============================================================================
cat("Testing core R packages...\n")
cat("----------------------------------------------------------\n")

core_packages <- c(
  "Seurat",
  "Signac", 
  "Matrix",
  "dplyr",
  "ggplot2",
  "glmnet",
  "ranger",
  "RANN"
)

for (pkg in core_packages) {
  if (test_package(pkg, required = TRUE)) {
    tests_passed <- tests_passed + 1
  } else {
    tests_failed <- tests_failed + 1
  }
}

# ============================================================================
# Test Bioconductor Packages
# ============================================================================
cat("\nTesting Bioconductor packages...\n")
cat("----------------------------------------------------------\n")

bioc_packages <- c(
  "GenomicRanges",
  "rtracklayer",
  "BiocGenerics"
)

for (pkg in bioc_packages) {
  if (test_package(pkg, required = TRUE)) {
    tests_passed <- tests_passed + 1
  } else {
    tests_failed <- tests_failed + 1
  }
}

# ============================================================================
# Test Genome Annotation Packages
# ============================================================================
cat("\nTesting genome annotation packages (species-specific)...\n")
cat("----------------------------------------------------------\n")

# Mouse
if (test_package("EnsDb.Mmusculus.v79", required = FALSE)) {
  tests_passed <- tests_passed + 1
  cat("    → Mouse (mm10) datasets supported\n")
} else {
  warnings_list <- c(warnings_list, "EnsDb.Mmusculus.v79 not installed - mouse data won't work")
}

# Human
if (test_package("EnsDb.Hsapiens.v86", required = FALSE)) {
  tests_passed <- tests_passed + 1
  cat("    → Human (hg38) datasets supported\n")
} else {
  warnings_list <- c(warnings_list, "EnsDb.Hsapiens.v86 not installed - human data won't work")
}

# ============================================================================
# Test Additional R Packages
# ============================================================================
cat("\nTesting additional R packages...\n")
cat("----------------------------------------------------------\n")

additional_packages <- c("patchwork", "viridis", "caret", "reticulate")
for (pkg in additional_packages) {
  if (test_package(pkg, required = TRUE)) {
    tests_passed <- tests_passed + 1
  } else {
    tests_failed <- tests_failed + 1
  }
}

# ============================================================================
# Test Python Environment (Required for neural networks)
# ============================================================================
cat("\nTesting Python environment (required for neural networks & autoencoders)...\n")
cat("----------------------------------------------------------\n")

python_available <- FALSE
tensorflow_available <- FALSE

if (requireNamespace("reticulate", quietly = TRUE)) {
  cat("  [✓] reticulate package found\n")
  tests_passed <- tests_passed + 1
  
  # Try to check Python
  tryCatch({
    py_config <- reticulate::py_config()
    cat(sprintf("  [✓] Python: %s\n", py_config$python))
    python_available <- TRUE
    
    # Check for TensorFlow (required for neural networks)
    tf_available <- tryCatch({
      reticulate::py_module_available("tensorflow")
    }, error = function(e) FALSE)
    
    if (tf_available) {
      cat("  [✓] TensorFlow available (neural networks supported)\n")
      tensorflow_available <- TRUE
      tests_passed <- tests_passed + 1
    } else {
      cat("  [!] TensorFlow NOT available - neural network training won't work\n")
      warnings_list <- c(warnings_list, "TensorFlow not installed - run: pip install tensorflow")
    }
    
    # Check for scvi-tools (optional for autoencoders)
    scvi_available <- tryCatch({
      reticulate::py_module_available("scvi")
    }, error = function(e) FALSE)
    
    if (scvi_available) {
      cat("  [✓] scvi-tools available (autoencoder methods supported)\n")
    } else {
      cat("  [!] scvi-tools NOT available - autoencoder methods won't work\n")
      warnings_list <- c(warnings_list, "scvi-tools not installed - scVI/PeakVI/MultiVI methods won't be available")
    }
    
  }, error = function(e) {
    cat(sprintf("  [!] Python not configured: %s\n", e$message))
    warnings_list <- c(warnings_list, "Python not configured - neural networks and autoencoders won't be available")
  })
} else {
  cat("  [✗] reticulate package NOT FOUND\n")
  tests_failed <- tests_failed + 1
  warnings_list <- c(warnings_list, "reticulate not installed - Python integration won't work")
}

# ============================================================================
# Test Configuration Files
# ============================================================================
cat("\nTesting configuration templates...\n")
cat("----------------------------------------------------------\n")

config_templates <- c(
  "combined/config_template.R",
  "rna_only/config_template.R",
  "atac_only/config_template.R"
)

for (template in config_templates) {
  if (file.exists(template)) {
    # Try to source the template
    result <- tryCatch({
      source(template, local = TRUE)
      TRUE
    }, error = function(e) {
      cat(sprintf("  [✗] %s - Error: %s\n", template, e$message))
      FALSE
    })
    
    if (result) {
      cat(sprintf("  [✓] %s\n", template))
      tests_passed <- tests_passed + 1
    } else {
      tests_failed <- tests_failed + 1
    }
  } else {
    cat(sprintf("  [✗] %s - NOT FOUND\n", template))
    tests_failed <- tests_failed + 1
  }
}

# ============================================================================
# Summary
# ============================================================================
cat("\n")
cat("============================================================\n")
cat("                    Test Summary                            \n")
cat("============================================================\n")
cat(sprintf("  Tests passed:  %d\n", tests_passed))
cat(sprintf("  Tests failed:  %d\n", tests_failed))
cat(sprintf("  Warnings:      %d\n", length(warnings_list)))

if (length(warnings_list) > 0) {
  cat("\nWarnings:\n")
  for (w in warnings_list) {
    cat(sprintf("  • %s\n", w))
  }
}

cat("\n")
if (tests_failed == 0) {
  cat("✓ Core installation successful! You can run the pipeline.\n")
  if (length(warnings_list) > 0) {
    cat("  Note: Some optional features are not available (see warnings above).\n")
  }
} else {
  cat("✗ Installation incomplete. Please install missing required packages:\n")
  cat("  source('combined/install_packages.R')\n")
}

cat("\n")
cat("Next steps:\n")
cat("  1. Copy config_template.R to config.R in your chosen pipeline directory\n")
cat("  2. Edit config.R with your dataset paths\n")
cat("  3. Run: Rscript run_pipeline.R\n")
cat("\n")
