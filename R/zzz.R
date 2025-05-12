.onAttach <- function(libname, pkgname) {
  deps <- c(
    "AnnotationDbi", "BiocParallel", "GENIE3", "INetTool", "Matrix",
    "RColorBrewer", "ReactomePA", "STRINGdb", "Seurat",
    "SingleCellExperiment", "SummarizedExperiment", "clusterProfiler",
    "distributions3", "doParallel", "dplyr", "fmsb", "ggraph", "ggplot2",
    "gridExtra", "igraph", "jsonlite", "org.Hs.eg.db", "parallel",
    "patchwork", "pROC", "rentrez", "reticulate", "robin", "scales", "tidyr"
  )

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  options(repos = BiocManager::repositories())

  for (pkg in deps) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      packageStartupMessage("scGraphVerse: installing: ", pkg)
      if (pkg %in% BiocManager::available()) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg)
      }
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, quietly = TRUE)
    )
  }

  if (requireNamespace("reticulate", quietly = TRUE)) {
    if (!reticulate::py_module_available("arboreto")) {
      packageStartupMessage("scGraphVerse: installing Python package 'arboreto'")
      reticulate::py_install("arboreto", pip = TRUE)
    }
  }

  if (!requireNamespace("ZILGM", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    packageStartupMessage("scGraphVerse: installing GitHub package 'bbeomjin/ZILGM'")
    remotes::install_github("bbeomjin/ZILGM", upgrade = "never")
  }

  if (!requireNamespace("learn2count", quietly = TRUE)) {
    packageStartupMessage("scGraphVerse: installing PCzinb ('drisso/learn2count')")
    BiocManager::install("drisso/learn2count", ask = FALSE, update = FALSE)
  }

  if (!requireNamespace("JRF", quietly = TRUE)) {
    packageStartupMessage("scGraphVerse: installing archived CRAN package 'JRF'")
    install.packages(
      "https://cran.r-project.org/src/contrib/Archive/JRF/JRF.1-4.tar.gz",
      repos = NULL,
      type  = "source"
    )
  }
  suppressPackageStartupMessages(
    library(JRF, character.only = TRUE, quietly = TRUE)
  )

  msg <- paste(
    cli::rule(center = "✨ Welcome to scGraphVerse ✨", line = 2),
    "",
    "🚀 Network analysis tools for single-cell data.",
    "",
    "This work is supported by the project:",
    "  National Centre for HPC, Big Data and Quantum Computing",
    "  Funded by European Union – Next Generation EU – CN00000013 🇪🇺",
    "  CUP: B93C22000620006",
    "",
    cli::rule(line = 2),
    sep = "\n"
  )
  packageStartupMessage(msg)
}

.onLoad <- function(libname, pkgname) {
}

.onUnload <- function(libpath) {
  library.dynam.unload("scGraphVerse", libpath)
}

