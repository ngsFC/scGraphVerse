# scGraphVerse

<div align="center">
  <img src="./man/figures/logo.png" alt="scGraphVerse logo" width="200"/>
  <h3>Comprehensive Gene Regulatory Network Analysis for Single-Cell Data</h3>
  
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
  [![R-CMD-check](https://github.com/ngsFC/scGraphVerse/workflows/R-CMD-check/badge.svg)](https://github.com/ngsFC/scGraphVerse/actions)
  [![Codecov test coverage](https://codecov.io/gh/ngsFC/scGraphVerse/branch/main/graph/badge.svg)](https://codecov.io/gh/ngsFC/scGraphVerse?branch=main)
  [![Bioc release status](http://www.bioconductor.org/shields/build/release/bioc/scGraphVerse.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/scGraphVerse)
  
  **Methods:** GENIE3 | GRNBoost2 | ZILGM | JRF | PCzinb
</div>

---

## Overview

**scGraphVerse** is a comprehensive R/Bioconductor package for inferring, evaluating, and visualizing gene regulatory networks (GRNs) from single-cell RNA sequencing data. It provides an integrated framework with multiple inference algorithms, consensus construction, and rich visualizations — all optimized for single-cell expression analysis.

### Key Features

- **🔬 Multiple Inference Methods**: GENIE3, GRNBoost2, ZILGM, JRF, PCzinb
- **🤝 Consensus Networks**: Voting, union, INet, and weighted consensus
- **📊 Comprehensive Evaluation**: ROC curves, AUC, F1-score, community analysis
- **🎨 Rich Visualizations**: Interactive networks, performance plots, publication-ready figures
- **🧬 Single-Cell Optimized**: Handles sparsity, batch effects, and zero-inflation
- **🔧 Flexible Integration**: Works with SingleCellExperiment, Seurat, and matrix objects

## Installation

### Bioconductor (Recommended)

```r
# Install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scGraphVerse")
```

### Development Version

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
- All R dependencies installed automatically

## Quick Start

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

# Detect communities and pathways
communities <- community_path(consensus)
```

## Method Overview

### Inference Algorithms

| Method | Description | Best For |
|--------|-------------|----------|
| **GENIE3** | Tree-based ensemble learning | General-purpose, robust |
| **GRNBoost2** | Gradient boosting with Dask | Large-scale datasets |
| **ZILGM** | Zero-inflated Gaussian graphical models | Sparse, zero-inflated data |
| **JRF** | Joint Random Forests | Multi-condition analysis |
| **PCzinb** | Partial correlation with ZINB | Direct regulatory relationships |

### Consensus Methods

- **Voting**: Democratic edge selection (majority vote)
- **Union**: Comprehensive edge inclusion (any method)
- **INet**: Intelligent network integration using INet-Tool
- **Weighted**: Score-based combination of edge weights

## Workflow

1. **Data Preparation**: Load single-cell RNA-seq count matrices
2. **Network Inference**: Apply multiple algorithms in parallel
3. **Consensus Building**: Integrate results using voting, union, or weighted methods
4. **Evaluation**: Assess performance using ROC curves, AUC, and network topology
5. **Visualization**: Generate publication-ready network plots and community analysis
6. **Biological Interpretation**: Pathway enrichment and functional analysis

## Examples

### Basic Network Inference

```r
# Single method inference
networks <- infer_networks(
  count_matrices_list = count_matrices,
  method = "GENIE3",
  nCores = 4
)

# Multiple methods
methods <- c("GENIE3", "ZILGM", "PCzinb")
all_networks <- lapply(methods, function(m) {
  infer_networks(count_matrices_list = count_matrices, method = m)
})
```

### Consensus Analysis

```r
# Create different consensus types
consensus_vote <- create_consensus(binary_adj, method = "vote")
consensus_union <- create_consensus(binary_adj, method = "union")
consensus_inet <- create_consensus(binary_adj, method = "INet")

# Compare consensus methods
comparison <- compare_consensus(
  consensus_matrix = consensus_vote,
  reference_matrix = adj_truth
)
```

### Performance Evaluation

```r
# ROC curve analysis
roc_results <- plotROC(
  swadj_list,
  adj_truth,
  plot_title = "Network Performance Comparison"
)

# Precision scores
precision_scores <- pscores(adj_truth, binary_adj)

# Community similarity
comm_similarity <- community_similarity(
  reference_communities,
  list(inferred_communities)
)
```

### Simulation and Benchmarking

```r
# Generate synthetic data
sim_data <- zinb_simdata(
  n = 50,           # cells per batch
  p = 100,          # genes
  B = adj_truth,    # ground truth network
  mu_range = list(c(1, 4), c(1, 7)),
  kmat = 2          # number of batches
)

# Benchmark methods
benchmark_results <- benchmark_methods(
  sim_data,
  methods = c("GENIE3", "ZILGM", "JRF"),
  ground_truth = adj_truth
)
```

## Documentation

- **Package Website**: https://ngsfc.github.io/scGraphVerse/
- **Vignettes**: 
  - [Quick Start Guide](https://ngsfc.github.io/scGraphVerse/articles/quickstart.html)
  - [Simulation Study](https://ngsfc.github.io/scGraphVerse/articles/simulation_study.html)
  - [Multi-Method Comparison](https://ngsfc.github.io/scGraphVerse/articles/case_study.html)
- **Reference Manual**: https://ngsfc.github.io/scGraphVerse/reference/

## Citation

```r
citation("scGraphVerse")
```

Please also cite the original methods:

1. **GENIE3**: Huynh-Thu et al. (2010). *PLOS ONE* 5(9):e12776
2. **GRNBoost2**: Moerman et al. (2019). *Bioinformatics* 35(12):2159-61
3. **ZILGM**: Park et al. (2021). *Statistical Analysis and Data Mining* 37(18):3085-3092
4. **JRF**: Petralia et al. (2015). *Journal of Proteome Research* 31(12):i197-i205
5. **PCzinb**: Nguyen et al. (2023). *Ann. Appl. Stat.* 17(3):2555-73
6. **INet-Tool**: Policastro et al. (2025). *Comput Stat* 40, 1517–1539

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

- **Bug Reports**: [GitHub Issues](https://github.com/ngsFC/scGraphVerse/issues)
- **Feature Requests**: [GitHub Issues](https://github.com/ngsFC/scGraphVerse/issues)
- **Pull Requests**: [GitHub Pull Requests](https://github.com/ngsFC/scGraphVerse/pulls)

## License

**scGraphVerse** is licensed under **GPL (≥ 2)**.

### Integrated Code Attribution

This package includes adapted functions from:
- **ZILGM** (Park et al., 2021) - GPL-2 license
- **JRF** (Petralia et al., 2015) - GPL (≥ 2) license

All integrated code maintains proper attribution and copyright notices.

## Funding

This work is supported by the **National Centre for HPC, Big Data and Quantum Computing**
- **Funded by**: European Union – Next Generation EU – CN00000013
- **CUP**: B93C22000620006

---

<div align="center">
  <strong>Happy Network Inference! 🧬📊</strong>
</div>