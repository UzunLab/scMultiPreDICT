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
#   source("install_packages.R")                  # From R console
#   Rscript install_packages.R                    # Core packages only
#   Rscript install_packages.R --species=mouse    # Add mouse annotations
#   Rscript install_packages.R --species=human    # Add human annotations
#   Rscript install_packages.R --species=skip     # Explicitly skip species packages
#   Rscript install_packages.R --lib-path-mode=auto|conda|default
#
# Recommended clean conda environment:
#   conda create -n scmulti -y -c conda-forge -c bioconda --strict-channel-priority \
#     python=3.10 r-base=4.3 sed coreutils grep gawk findutils zlib libzlib libxml2 libcurl openssl git
#   conda activate scmulti
#
# ============================================================================

cat("============================================================\n")
cat("      scMultiPreDICT R Package Installation                 \n")
cat("============================================================\n\n")

CRAN_REPO <- "https://cloud.r-project.org"
options(repos = c(CRAN = CRAN_REPO))
r_version <- getRversion()

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
species_arg <- NULL
lib_path_mode <- "auto"
for (arg in args) {
  if (grepl("^--species=", arg)) {
    species_arg <- tolower(sub("^--species=", "", arg))
  } else if (grepl("^--lib-path-mode=", arg)) {
    lib_path_mode <- tolower(sub("^--lib-path-mode=", "", arg))
  }
}
if (!is.null(species_arg) && !(species_arg %in% c("mouse", "human", "skip"))) {
  cat(sprintf("WARNING: Unknown --species value '%s'. Using 'skip'.\n", species_arg))
  species_arg <- "skip"
}
if (!(lib_path_mode %in% c("auto", "conda", "default"))) {
  cat(sprintf("WARNING: Unknown --lib-path-mode value '%s'. Using 'auto'.\n", lib_path_mode))
  lib_path_mode <- "auto"
}

cat(sprintf("R version detected: %s\n", as.character(r_version)))
cat(sprintf("Platform: %s\n", R.version$platform))

conda_bin <- Sys.which("conda")
conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
conda_detected <- nzchar(conda_bin)
conda_available <- FALSE
if (conda_detected) {
  cat(sprintf("Conda detected: %s\n", conda_bin))
  if (nzchar(conda_prefix)) {
    cat(sprintf("Conda target environment: %s\n", conda_prefix))
    conda_history <- file.path(conda_prefix, "conda-meta", "history")
    if (file.exists(conda_history)) {
      conda_available <- TRUE
    } else {
      cat("WARNING: Active conda prefix is missing 'conda-meta/history'.\n")
      cat("         Conda fallback installs are disabled for this run.\n")
      cat("         Recreate this environment before relying on conda fallback.\n")
      cat("         Example:\n")
      cat("           conda create -n scmulti -y -c conda-forge -c bioconda --strict-channel-priority \\\n")
      cat("             python=3.10 r-base=4.3 sed coreutils grep gawk findutils zlib libzlib libxml2 libcurl openssl git\n")
      cat("           conda activate scmulti\n")
    }
  } else {
    cat("WARNING: No active conda prefix (CONDA_PREFIX is empty).\n")
    cat("         Conda fallback installs are disabled to avoid polluting base env.\n")
  }
} else {
  cat("Conda not found in PATH. Conda fallback will be skipped.\n")
}

configure_library_paths <- function(conda_prefix_path, mode = "auto") {
  current_paths <- .libPaths()

  if (mode == "default") {
    cat("Library path mode: default (keeping existing .libPaths()).\n")
    cat(sprintf("Active library paths:\n  %s\n", paste(current_paths, collapse = "\n  ")))
    return(invisible(FALSE))
  }

  if (!nzchar(conda_prefix_path)) {
    cat("No active conda prefix detected. Keeping current .libPaths().\n")
    cat(sprintf("Active library paths:\n  %s\n", paste(current_paths, collapse = "\n  ")))
    return(invisible(FALSE))
  }

  conda_r_lib <- file.path(conda_prefix_path, "lib", "R", "library")
  if (!dir.exists(conda_r_lib)) {
    cat("Conda prefix detected, but conda R library path is missing.\n")
    cat(sprintf("Expected path: %s\n", conda_r_lib))
    cat("Keeping current .libPaths().\n")
    cat(sprintf("Active library paths:\n  %s\n", paste(current_paths, collapse = "\n  ")))
    return(invisible(FALSE))
  }

  conda_root_norm <- normalizePath(conda_prefix_path, mustWork = FALSE)
  current_norm <- vapply(current_paths, normalizePath, FUN.VALUE = character(1), mustWork = FALSE)
  has_non_conda_paths <- any(!startsWith(current_norm, conda_root_norm))

  should_isolate <- identical(mode, "conda") || (identical(mode, "auto") && has_non_conda_paths)

  if (!should_isolate) {
    cat(sprintf("Library path mode: %s (no isolation needed).\n", mode))
    cat(sprintf("Active library paths:\n  %s\n", paste(current_paths, collapse = "\n  ")))
    return(invisible(FALSE))
  }

  # In conda environments, avoid mixing user-level compiled packages from
  # ~/R/... with conda binaries; this is a common source of .so load errors.
  .libPaths(conda_r_lib)
  Sys.setenv(R_LIBS_USER = conda_r_lib)

  cat(sprintf("Library path mode: %s (configured conda-isolated .libPaths()).\n", mode))
  cat(sprintf("Previous library paths:\n  %s\n", paste(current_paths, collapse = "\n  ")))
  cat(sprintf("Current library paths:\n  %s\n", paste(.libPaths(), collapse = "\n  ")))
  invisible(TRUE)
}

configure_library_paths(conda_prefix, mode = lib_path_mode)
cat("\n")

# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------

pkg_status <- function(pkg) {
  installed <- pkg %in% rownames(installed.packages())
  if (!installed) {
    return(list(
      installed = FALSE,
      loadable = FALSE,
      version = NA_character_,
      error = NULL
    ))
  }

  version <- tryCatch(
    as.character(packageVersion(pkg)),
    error = function(e) NA_character_
  )

  load_error <- NULL
  loadable <- TRUE
  tryCatch(
    suppressPackageStartupMessages(loadNamespace(pkg)),
    error = function(e) {
      loadable <<- FALSE
      load_error <<- conditionMessage(e)
    }
  )

  list(
    installed = TRUE,
    loadable = loadable,
    version = version,
    error = load_error
  )
}

ensure_remotes <- function() {
  st <- pkg_status("remotes")
  if (st$loadable) {
    return(TRUE)
  }
  cat("    Installing helper package: remotes\n")
  tryCatch(
    install.packages("remotes", repos = CRAN_REPO, quiet = FALSE),
    error = function(e) cat(sprintf("    remotes install error: %s\n", e$message))
  )

  st <- pkg_status("remotes")
  if (!st$loadable && conda_available) {
    cat("    CRAN remotes install unavailable; trying conda fallback.\n")
    conda_status <- tryCatch(
      system2(conda_bin, args = c("install", "-y", "-p", conda_prefix, "conda-forge::r-remotes")),
      error = function(e) {
        cat(sprintf("    conda remotes install error: %s\n", e$message))
        1L
      }
    )
    if (!identical(conda_status, 0L)) {
      cat(sprintf("    conda remotes install exited with status %s.\n", conda_status))
    }
  }
  pkg_status("remotes")$loadable
}

ensure_biocmanager <- function() {
  st <- pkg_status("BiocManager")
  if (st$loadable) {
    return(TRUE)
  }
  cat("    Installing helper package: BiocManager\n")
  tryCatch(
    install.packages("BiocManager", repos = CRAN_REPO, quiet = FALSE),
    error = function(e) cat(sprintf("    BiocManager install error: %s\n", e$message))
  )

  st <- pkg_status("BiocManager")
  if (!st$loadable && conda_available) {
    cat("    CRAN BiocManager install unavailable; trying conda fallback.\n")
    conda_status <- tryCatch(
      system2(conda_bin, args = c("install", "-y", "-p", conda_prefix, "conda-forge::r-biocmanager")),
      error = function(e) {
        cat(sprintf("    conda BiocManager install error: %s\n", e$message))
        1L
      }
    )
    if (!identical(conda_status, 0L)) {
      cat(sprintf("    conda BiocManager install exited with status %s.\n", conda_status))
    }
  }
  pkg_status("BiocManager")$loadable
}

install_from_cran <- function(pkg) {
  install.packages(pkg, repos = CRAN_REPO, quiet = FALSE)
}

install_from_cran_archive <- function(pkg, version) {
  if (!ensure_remotes()) {
    cat("    Cannot install CRAN archive package because remotes is unavailable.\n")
    return(invisible(FALSE))
  }
  remotes::install_version(
    package = pkg,
    version = version,
    repos = CRAN_REPO,
    upgrade = "never",
    quiet = FALSE
  )
  invisible(TRUE)
}

install_from_github <- function(repo, ref = NULL) {
  if (!ensure_remotes()) {
    cat("    Cannot install from GitHub because remotes is unavailable.\n")
    return(invisible(FALSE))
  }
  if (is.null(ref)) {
    remotes::install_github(repo, quiet = FALSE, upgrade = "never")
  } else {
    remotes::install_github(repo, ref = ref, quiet = FALSE, upgrade = "never")
  }
  invisible(TRUE)
}

install_from_bioc <- function(pkg, force = FALSE) {
  if (!ensure_biocmanager()) {
    cat("    Cannot install via BiocManager because BiocManager is unavailable.\n")
    return(invisible(FALSE))
  }
  BiocManager::install(pkg, update = FALSE, ask = FALSE, force = force)
  invisible(TRUE)
}

install_from_conda <- function(spec) {
  if (!conda_available) {
    cat("    Conda not available. Skipping conda fallback.\n")
    return(invisible(FALSE))
  }

  conda_args <- c("install", "-y")
  if (nzchar(conda_prefix)) {
    conda_args <- c(conda_args, "-p", conda_prefix)
  }
  conda_args <- c(conda_args, spec)

  status <- tryCatch(
    system2(conda_bin, args = conda_args),
    error = function(e) {
      cat(sprintf("    Conda install error for '%s': %s\n", spec, e$message))
      1L
    }
  )

  if (!identical(status, 0L)) {
    cat(sprintf("    Conda command exited with status %s for '%s'.\n", status, spec))
  }
  invisible(identical(status, 0L))
}

make_method <- function(label, fn) {
  list(label = label, fn = fn)
}

append_method <- function(methods, method) {
  methods[[length(methods) + 1]] <- method
  methods
}

install_with_methods <- function(pkg, methods) {
  st <- pkg_status(pkg)
  if (st$loadable) {
    cat(sprintf("  Already installed: %s (%s)\n", pkg, st$version))
    return(TRUE)
  }

  if (st$installed && !st$loadable) {
    cat(sprintf("  %s is installed but not loadable.\n", pkg))
    cat(sprintf("    Load error: %s\n", st$error))
  }

  for (method in methods) {
    cat(sprintf("  Installing: %s via %s\n", pkg, method$label))
    tryCatch(
      method$fn(),
      error = function(e) cat(sprintf("    %s error: %s\n", method$label, e$message))
    )

    st <- pkg_status(pkg)
    if (st$loadable) {
      cat(sprintf("    Success: %s (%s)\n", pkg, st$version))
      return(TRUE)
    }

    if (st$installed && !st$loadable) {
      cat(sprintf("    Installed but still not loadable: %s\n", st$error))
    } else {
      cat(sprintf("    %s still not installed.\n", pkg))
    }
  }

  cat(sprintf("  Failed: %s\n", pkg))
  FALSE
}

# ----------------------------------------------------------------------------
# Package groups and fallback maps
# ----------------------------------------------------------------------------

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

  # Machine learning and utilities
  "glmnet",
  "ranger",
  "RANN",
  "irlba",
  "reticulate",
  "png",
  "parallel",
  "doParallel",
  "foreach",
  "matrixStats"
)

bioc_packages <- c(
  "GenomicRanges",
  "GenomeInfoDb",
  "IRanges",
  "S4Vectors",
  "rtracklayer"
)

conda_map <- c(
  # CRAN
  dplyr = "conda-forge::r-dplyr",
  tidyr = "conda-forge::r-tidyr",
  readr = "conda-forge::r-readr",
  purrr = "conda-forge::r-purrr",
  stringr = "conda-forge::r-stringr",
  tibble = "conda-forge::r-tibble",
  forcats = "conda-forge::r-forcats",
  ggplot2 = "conda-forge::r-ggplot2",
  patchwork = "conda-forge::r-patchwork",
  viridis = "conda-forge::r-viridis",
  scales = "conda-forge::r-scales",
  RColorBrewer = "conda-forge::r-rcolorbrewer",
  cowplot = "conda-forge::r-cowplot",
  RANN = "conda-forge::r-rann",
  irlba = "conda-forge::r-irlba",
  doParallel = "conda-forge::r-doparallel",
  foreach = "conda-forge::r-foreach",
  matrixStats = "conda-forge::r-matrixstats",
  Matrix = "conda-forge::r-matrix",
  MASS = "conda-forge::r-mass",
  caret = "conda-forge::r-caret",
  Seurat = "conda-forge::r-seurat",
  Signac = "bioconda::r-signac",
  png = "conda-forge::r-png",
  ggrepel = "conda-forge::r-ggrepel",
  reshape2 = "conda-forge::r-reshape2",
  glmnet = "conda-forge::r-glmnet",
  ranger = "conda-forge::r-ranger",
  reticulate = "conda-forge::r-reticulate",
  # Bioconductor
  GenomicRanges = "bioconda::bioconductor-genomicranges",
  GenomeInfoDb = "bioconda::bioconductor-genomeinfodb",
  IRanges = "bioconda::bioconductor-iranges",
  S4Vectors = "bioconda::bioconductor-s4vectors",
  rtracklayer = "bioconda::bioconductor-rtracklayer"
)

get_conda_spec <- function(pkg) {
  spec <- unname(conda_map[pkg])
  if (length(spec) == 0L || is.na(spec)) {
    return(NULL)
  }
  spec
}

# ----------------------------------------------------------------------------
# Installation steps
# ----------------------------------------------------------------------------

failed_packages <- c()

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  methods <- list()

  if (pkg == "ggrepel" && r_version < "4.5.0") {
    methods <- append_method(
      methods,
      make_method("CRAN archive (ggrepel 0.9.6)", function() {
        install_from_cran_archive("ggrepel", "0.9.6")
      })
    )
  }

  if (pkg == "reshape2") {
    methods <- append_method(
      methods,
      make_method("CRAN archive (reshape2 1.4.4)", function() {
        install_from_cran_archive("reshape2", "1.4.4")
      })
    )
  }

  methods <- append_method(
    methods,
    make_method("CRAN", function() install_from_cran(pkg))
  )

  conda_spec <- get_conda_spec(pkg)
  if (!is.null(conda_spec)) {
    methods <- append_method(
      methods,
      make_method(sprintf("conda (%s)", conda_spec), function() install_from_conda(conda_spec))
    )
  }

  ok <- install_with_methods(pkg, methods)
  if (!ok) {
    failed_packages <- unique(c(failed_packages, pkg))
  }
}

cat("\nInstalling Seurat stack...\n")
for (pkg in c("Matrix", "MASS", "caret")) {
  methods <- list(
    make_method("CRAN", function() install_from_cran(pkg))
  )
  conda_spec <- get_conda_spec(pkg)
  if (!is.null(conda_spec)) {
    methods <- append_method(
      methods,
      make_method(sprintf("conda (%s)", conda_spec), function() install_from_conda(conda_spec))
    )
  }
  ok <- install_with_methods(pkg, methods)
  if (!ok) {
    failed_packages <- unique(c(failed_packages, pkg))
  }
}

seurat_methods <- list(
  make_method("CRAN", function() install_from_cran("Seurat"))
)
seurat_conda <- get_conda_spec("Seurat")
if (!is.null(seurat_conda)) {
  seurat_methods <- append_method(
    seurat_methods,
    make_method(sprintf("conda (%s)", seurat_conda), function() install_from_conda(seurat_conda))
  )
}
seurat_methods <- append_method(
  seurat_methods,
  make_method("GitHub (satijalab/seurat@seurat5)", function() {
    install_from_github("satijalab/seurat", ref = "seurat5")
  })
)

if (!install_with_methods("Seurat", seurat_methods)) {
  failed_packages <- unique(c(failed_packages, "Seurat"))
}

cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
  methods <- list()

  conda_spec <- get_conda_spec(pkg)
  if (!is.null(conda_spec)) {
    methods <- append_method(
      methods,
      make_method(sprintf("conda (%s)", conda_spec), function() install_from_conda(conda_spec))
    )
  }

  methods <- append_method(
    methods,
    make_method("BiocManager", function() install_from_bioc(pkg, force = FALSE))
  )
  methods <- append_method(
    methods,
    make_method("BiocManager (force reinstall)", function() install_from_bioc(pkg, force = TRUE))
  )

  ok <- install_with_methods(pkg, methods)
  if (!ok) {
    failed_packages <- unique(c(failed_packages, pkg))
  }
}

signac_methods <- list(
  make_method("CRAN", function() install_from_cran("Signac"))
)
signac_conda <- get_conda_spec("Signac")
if (!is.null(signac_conda)) {
  signac_methods <- append_method(
    signac_methods,
    make_method(sprintf("conda (%s)", signac_conda), function() install_from_conda(signac_conda))
  )
}
signac_methods <- append_method(
  signac_methods,
  make_method("GitHub (stuart-lab/signac)", function() install_from_github("stuart-lab/signac"))
)

cat("\nInstalling Signac...\n")
if (!install_with_methods("Signac", signac_methods)) {
  failed_packages <- unique(c(failed_packages, "Signac"))
}

# ----------------------------------------------------------------------------
# Species-specific annotation packages
# ----------------------------------------------------------------------------

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

if (!is.null(species_arg)) {
  species_choice <- species_arg
  cat(sprintf("Species specified via command line: %s\n", species_choice))
} else if (interactive()) {
  species_choice <- tolower(readline(prompt = "Install annotation packages? (mouse/human/skip): "))
} else {
  species_choice <- "skip"
  cat("No --species argument provided. Skipping species-specific packages.\n")
}

if (species_choice == "mouse") {
  cat("Installing mouse annotation packages...\n")
  tryCatch(
    {
      install_from_bioc("EnsDb.Mmusculus.v79", force = FALSE)
      install_from_bioc("BSgenome.Mmusculus.UCSC.mm10", force = FALSE)
    },
    error = function(e) cat(sprintf("  Mouse package install error: %s\n", e$message))
  )
} else if (species_choice == "human") {
  cat("Installing human annotation packages...\n")
  tryCatch(
    {
      install_from_bioc("EnsDb.Hsapiens.v86", force = FALSE)
      install_from_bioc("BSgenome.Hsapiens.UCSC.hg38", force = FALSE)
    },
    error = function(e) cat(sprintf("  Human package install error: %s\n", e$message))
  )
} else {
  cat("Skipping species-specific packages. Install manually when needed.\n")
}

# ----------------------------------------------------------------------------
# Installation verification
# ----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("Installation Verification\n")
cat("============================================================\n")

all_required_packages <- unique(c(
  cran_packages,
  "Matrix", "MASS", "caret", "Seurat", "Signac",
  bioc_packages
))

not_installed <- c()
not_loadable <- c()
load_errors <- list()

for (pkg in all_required_packages) {
  st <- pkg_status(pkg)
  if (!st$installed) {
    not_installed <- c(not_installed, pkg)
    next
  }
  if (!st$loadable) {
    not_loadable <- c(not_loadable, pkg)
    load_errors[[pkg]] <- st$error
  }
}

if (length(not_installed) == 0 && length(not_loadable) == 0) {
  cat("\nAll core packages are installed and loadable.\n")
} else {
  cat("\nWARNING: Some packages are still unavailable.\n")

  if (length(not_installed) > 0) {
    cat("\nNot installed:\n")
    cat(paste("  -", unique(not_installed), collapse = "\n"), "\n")
  }

  if (length(not_loadable) > 0) {
    cat("\nInstalled but not loadable:\n")
    for (pkg in unique(not_loadable)) {
      cat(sprintf("  - %s\n", pkg))
      err <- load_errors[[pkg]]
      if (!is.null(err) && nzchar(err)) {
        cat(sprintf("      Load error: %s\n", err))
      }
    }
  }
}

if (length(failed_packages) > 0) {
  cat("\nPackages that failed all automated install methods:\n")
  cat(paste("  -", unique(failed_packages), collapse = "\n"), "\n")
}

# ----------------------------------------------------------------------------
# Session information
# ----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("Session Information\n")
cat("============================================================\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform: %s\n", R.version$platform))
cat(sprintf("Library paths:\n  %s\n", paste(.libPaths(), collapse = "\n  ")))

for (pkg in c("BiocManager", "Seurat", "Signac", "reticulate", "caret")) {
  st <- pkg_status(pkg)
  if (st$loadable) {
    cat(sprintf("%s version: %s\n", pkg, st$version))
  } else if (st$installed) {
    cat(sprintf("%s: INSTALLED but NOT LOADABLE\n", pkg))
  } else {
    cat(sprintf("%s: NOT INSTALLED\n", pkg))
  }
}

cat("\n============================================================\n")
cat("Installation Complete\n")
cat("============================================================\n")
