# scGraphVerse

<div align="center">
  <img src="./man/figures/logo.png" alt="scGraphVerse logo" width="200"/>
  <h3>Comprehensive Gene Regulatory Network Analysis for Single-Cell Data</h3>
  
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
  
</div>

---

## Overview

**scGraphVerse** is a comprehensive R package for inferring, evaluating, and visualizing gene regulatory networks (GRNs) from single-cell RNA sequencing data. It provides an integrated framework with multiple inference algorithms, consensus construction, and rich visualizations optimized for single-cell expression analysis.

### Key Features

- **🔬 Multiple Inference Methods**: GENIE3, GRNBoost2, ZILGM, JRF, PCzinb
- **🤝 Consensus Networks**: Voting, union, INet
- **📊 Comprehensive Evaluation**: ROC curves, AUC, F1-score, community analysis
- **🎨 Rich Visualizations**: Interactive networks, performance plots, publication-ready figures
- **🔧 Flexible Integration**: Works with SingleCellExperiment, Seurat, and matrix objects

## Installation

### Development Version

```r
# Install development version
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ngsFC/scGraphVerse")
```

### Inference Algorithms

| Method | Description |
|--------|-------------|
| **GENIE3** | Tree-based ensemble learning |
| **GRNBoost2** | Gradient boosting with Dask |
| **ZILGM** | Zero-inflated Gaussian graphical models | 
| **JRF** | Joint Random Forests |
| **PCzinb** | Partial correlation with ZINB |


## Documentation

- **Package Website**: https://ngsfc.github.io/scGraphVerse/
- **Vignettes**: 
  - [Simulation Study](https://ngsfc.github.io/scGraphVerse/articles/simulation_study.html)
  - [Case Study](https://ngsfc.github.io/scGraphVerse/articles/case_study.html)
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
