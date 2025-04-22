# scGraphVerse
               
<div align="center">
  <img src="logo.png" alt="Logo" width="200"/>
</div>

**scGraphVerse** is an R package for inferring, evaluating, and visualizing gene regulatory networks (GRNs) from single-cell RNA-seq data. It integrates multiple GRN inference methods, customizable thresholding strategies, performance metrics, and consensus network construction — all designed with flexibility, scalability, and visualization in mind.

---

## ✨ Features

- 🔍 Inference of GRNs from count matrices using methods like **GENIE3**, **GRNBoost2**, **ZILGM**, **PCzimb** and **JRF**.
- 🎯 Thresholding of weighted networks using shuffled matrix null models `cutoff_adjacency`.
- 🧠 Consensus network generation via `voting`, `union`, or the `INet` framework.
- 📊 Performance evaluation with ROC curves, classification metrics and topological measures.
- 🧩 Network visualization with `ggraph`.

---

<div align="center">
  <img src="gabstract.png" alt="Logo" width="700"/>
</div>

## 🧬 Installation

This package is not yet on CRAN or Bioconductor. You can install the development version directly from GitHub:

```r
# You need devtools to install from GitHub
install.packages("devtools")
devtools::install_github("ngsFC/scGraphVerse")
```

---
## Abstract BITS 2025

### Title
**scGraphVerse**: A Unified Framework for Network Inference, Evaluation, and Visualization from Single-Cell Gene Expression Data

### Authors
Francesco Cecere, Daniela De Canditiis, Annamaria Carissimo, Claudia Angelini

### Affiliations
Institute of Genetics and Biophysics (IGB) "Adriano Buzzati-Traverso", Consiglio Nazionale delle Ricerche (CNR), Naples, Italy
Istituto per le Applicazioni del Calcolo (IAC) "Mauro Picone", Consiglio Nazionale delle Ricerche (CNR), 80131 Naples, Italy

### Motivation
Gene regulatory network (GRN) inference from single-cell RNA sequencing (scRNA-seq) data is essential for understanding the molecular dynamics of a specific cellular type. However, the diversity of available network inference algorithms and methods and lack of standard evaluation pipelines, in particular for the integration of several experiments, pose challenges for benchmarking and reproducibility. Existing tools often focus on a single task (e.g., inference or visualization), are method-specific, or lack support for modular workflows across different data types and conditions. To address these limitations, we developed scGraphVerse, an R/Bioconductor package that provides a unified and extensible framework for GRN inference, evaluation, consensus building, and visualization from scRNA-seq data.

### Methods
scGraphVerse supports a modular workflow where users can input single or multiple gene expression matrices and apply a variety of inference algorithms. The package allows users to compare network structures generated from early, late, or joint data integration strategies. scGraphVerse includes utilities to compute network similarity, identify consensus edges across runs or methods, and evaluate networks using gold standard interactions. 
Taking common metrics such as TPR, FPR, Precision and others it can evaluate the performance assessment of the methods. scGraphVerse also facilitates downstream analysis such as community detection, pathway analysis and text mining approach to study gene-gene interaction in PubMed. Visualization modules provide publication-ready plots for network structure, ROC curves, and community similarity. All functions are compatible with standard Bioconductor data structures, ensuring interoperability with existing workflows.

### Results
We applied scGraphVerse to benchmark network inference methods across simulated and real single-cell datasets. In simulation, we evaluated performance across varying numbers of experiments (1 to 5 matrices), repeated runs (n = 10), and integration strategies (late, early, and joint) using GENIE3, GRNBoost2, and JRF. Using synthetic scRNA-seq data (n = 100 cells, k = 3 matrices), we observed that early integration with GENIE3 achieved the highest true positive rate (TPR ≈ 0.43) and F1-score (F1 ≈ 0.50) at smaller gene sets (p = 82), while late integration with GRNBoost2 maintained a low false positive rate (FPR ≈ 0.000–0.001) and high precision (Precision ≈ 0.92). However, JRF consistently balanced these trade-offs, achieving F1-scores around 0.30–0.33, MCC values near 0.30–0.40, and precision up to 0.90, with nearly zero FPR across all gene sizes. While JRF’s runtime increased with gene count (0.5 to 10 minutes), it provided strong accuracy and generalization, particularly for moderate-sized gene sets.

### References

GENIE3: Huynh-Thu et al., 2010. Inferring regulatory networks from expression data using tree-based methods. Bioinformatics.

GRNBoost2: Moerman et al., 2019. GRNBoost2 and Arboreto: Efficient and scalable inference of gene regulatory networks. Bioinformatics.

ZILGM: Zhang et al., 2021. A zero-inflated log-normal model for single-cell gene expression networks. Bioinformatics.

JRF: Shu et al., 2016. Joint estimation of multiple dependent Gaussian graphical models with applications to mouse genomics. Journal of Machine Learning Research.

### Abstract URL

### Supplemetary info

---

Progetto "National Centre for HPC, Big Data and Quantum Computing" founded by European Union - Next Generation EU - CN00000013 - CUP B93C22000620006

---
