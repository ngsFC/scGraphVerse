---
title: "scGraphVerse"
output: rmarkdown::github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

<div class="hero-section">
  <div class="container text-center">
    <img src="./man/figures/logo.png" alt="scGraphVerse logo" width="200" style="margin-bottom: 30px;"/>
    <h1>scGraphVerse</h1>
    <p class="lead">Comprehensive Gene Regulatory Network Analysis for Single-Cell Data</p>
    <div style="margin-top: 30px;">
      <span class="method-badge">GENIE3</span>
      <span class="method-badge">GRNBoost2</span>
      <span class="method-badge">ZILGM</span>
      <span class="method-badge">JRF</span>
      <span class="method-badge">PCzinb</span>
    </div>
  </div>
</div>

<div align="center">

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/ngsFC/scGraphVerse/workflows/R-CMD-check/badge.svg)](https://github.com/ngsFC/scGraphVerse/actions)
[![Codecov test coverage](https://codecov.io/gh/ngsFC/scGraphVerse/branch/main/graph/badge.svg)](https://codecov.io/gh/ngsFC/scGraphVerse?branch=main)
[![CRAN status](https://www.r-pkg.org/badges/version/scGraphVerse)](https://CRAN.R-project.org/package=scGraphVerse)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/scGraphVerse)](https://cran.r-project.org/package=scGraphVerse)
[![Bioc release status](http://www.bioconductor.org/shields/build/release/bioc/scGraphVerse.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/scGraphVerse)

</div>

**scGraphVerse** is a comprehensive R/Bioconductor package for inferring, evaluating, and visualizing gene regulatory networks (GRNs) from single-cell RNA sequencing data. It provides an integrated framework with multiple inference algorithms, enhanced parameter control, consensus construction, and rich visualizations — all optimized for single-cell expression analysis.

---

## ✨ Key Features

<div class="row">
  <div class="col-md-6">
    <div class="feature-box">
      <h3><i class="fas fa-network-wired text-primary"></i> Multiple Inference Methods</h3>
      <ul>
        <li><strong>GENIE3</strong> - Tree-based ensemble learning</li>
        <li><strong>GRNBoost2</strong> - Gradient boosting with Dask</li>
        <li><strong>ZILGM</strong> - Zero-inflated Gaussian graphical models</li>
        <li><strong>JRF</strong> - Joint Random Forests</li>
        <li><strong>PCzinb</strong> - Partial correlation with ZINB</li>
      </ul>
    </div>
  </div>
  <div class="col-md-6">
    <div class="feature-box">
      <h3><i class="fas fa-layer-group text-success"></i> Advanced Consensus</h3>
      <ul>
        <li><strong>Voting consensus</strong> - Democratic edge selection</li>
        <li><strong>Union consensus</strong> - Comprehensive edge inclusion</li>
        <li><strong>INet consensus</strong> - Intelligent network integration</li>
        <li><strong>Weighted consensus</strong> - Score-based combination</li>
      </ul>
    </div>
  </div>
</div>

<div class="row">
  <div class="col-md-6">
    <div class="feature-box">
      <h3><i class="fas fa-chart-bar text-info"></i> Comprehensive Evaluation</h3>
      <ul>
        <li><strong>Performance metrics</strong> - ROC, AUC, F1-score, MCC</li>
        <li><strong>Network topology</strong> - Modularity, clustering coefficients</li>
        <li><strong>Community detection</strong> - With pathway enrichment</li>
        <li><strong>Benchmarking</strong> - Against known networks</li>
      </ul>
    </div>
  </div>
  <div class="col-md-6">
    <div class="feature-box">
      <h3><i class="fas fa-paint-brush text-warning"></i> Rich Visualizations</h3>
      <ul>
        <li><strong>Interactive networks</strong> - Community highlighting</li>
        <li><strong>Performance plots</strong> - ROC curves, precision-recall</li>
        <li><strong>Consensus comparisons</strong> - Method benchmarking</li>
        <li><strong>Publication-ready</strong> - High-quality figures</li>
      </ul>
    </div>
  </div>
</div>

## 🔄 Workflow Overview

<div class="row text-center">
  <div class="col-md-2">
    <div class="workflow-step">
      <i class="fas fa-upload fa-2x text-primary"></i>
      <h5>1. Load Data</h5>
      <p>Single-cell RNA-seq matrices</p>
    </div>
  </div>
  <div class="col-md-1">
    <div class="workflow-arrow">→</div>
  </div>
  <div class="col-md-2">
    <div class="workflow-step">
      <i class="fas fa-network-wired fa-2x text-success"></i>
      <h5>2. Infer Networks</h5>
      <p>Multiple algorithms</p>
    </div>
  </div>
  <div class="col-md-1">
    <div class="workflow-arrow">→</div>
  </div>
  <div class="col-md-2">
    <div class="workflow-step">
      <i class="fas fa-layer-group fa-2x text-info"></i>
      <h5>3. Build Consensus</h5>
      <p>Integrate results</p>
    </div>
  </div>
  <div class="col-md-1">
    <div class="workflow-arrow">→</div>
  </div>
  <div class="col-md-2">
    <div class="workflow-step">
      <i class="fas fa-chart-line fa-2x text-warning"></i>
      <h5>4. Analyze & Visualize</h5>
      <p>Communities & pathways</p>
    </div>
  </div>
</div>

---

<div align="center">
  <img src="./man/figures/gabstract.png" alt="Graphical Abstract" width="700"/>
</div>

---

## 🚀 Quick Start

```r
library(scGraphVerse)

# Load example data
data("count_matrices")

# Infer networks using GENIE3
networks <- infer_networks(
  count_matrices_list = count_matrices,
  method = "GENIE3",
  nCores = 4
)

# Create consensus network
wadj <- generate_adjacency(networks)
swadj <- symmetrize(wadj, weight_function = "mean")
binary_adj <- cutoff_adjacency(
  count_matrices = count_matrices,
  weighted_adjm_list = swadj,
  quantile_threshold = 0.95
)
consensus <- create_consensus(binary_adj, method = "vote")

# Visualize network
plotg(list(consensus))

# Detect communities
communities <- community_path(consensus)
```

## 💾 Installation

### Bioconductor (Recommended)

```r
# Install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scGraphVerse")
```

### From GitHub (Development Version)

```r
# Install development version
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ngsFC/scGraphVerse")

# For GRNBoost2 support (optional)
library(scGraphVerse)
init_py(install_missing = TRUE)  # Installs Python dependencies
```

### System Requirements

- **R** ≥ 4.4.0
- **Python** ≥ 3.6 (optional, for GRNBoost2)
- **Dependencies**: All R dependencies installed automatically

---

---

## 📄 License & Attribution

**scGraphVerse** is licensed under **GPL (>= 2)** .

### Integrated Code Attribution

This package includes adapted functions from:

- **ZILGM** (Park et al., 2021) - GPL-2 license
- **JRF** (Petralia et al., 2015) - GPL (>= 2) license

All integrated code maintains proper attribution and copyright notices.

---

## 📚 Citation

```r
# Get citation information
citation("scGraphVerse")
```

Please also cite the original methods implemented:

1. **GENIE3**: Huynh-Thu et al. (2010). *PLOS ONE* 5(9):e12776
2. **GRNBoost2**: Moerman et al. (2019). *Bioinformatics* 35(12):2159-61
3. **ZILGM**: Park et al. (2021). *Statistical Analysis and Data Mining* 37(18):3085-3092
4. **JRF**: Petralia et al. (2015). *Journal of Proteome Research* 31(12):i197-i205
5. **PCzinb**: Nguyen et al. (2023). *Ann. Appl. Stat.* 17(3):2555-73
6. **INet-Tool**: Policastro et al (2025). *Comput Stat* 40, 1517–1539

---

---

## 💰 Funding

This work is supported by the **National Centre for HPC, Big Data and Quantum Computing** 
- **Funded by**: European Union – Next Generation EU – CN00000013 
- **CUP**: B93C22000620006

---

---

<div align="center">
  <p><strong>Happy Network Inference! 🧬📊</strong></p>
</div>
