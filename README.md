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

Progetto "National Centre for HPC, Big Data and Quantum Computing" founded by European Union - Next Generation EU - CN00000013 - CUP B93C22000620006

---
