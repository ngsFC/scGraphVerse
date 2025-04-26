# scGraphVerse News

## scGraphVerse 0.1.0 (Development Version)

### 🚀 Initial Release

- Released the first development version of `scGraphVerse`.
- Provides gene regulatory network (GRN) inference from single-cell RNA-seq data.
- Supports multiple network inference methods:
  - **GENIE3** (tree-based ensemble method)
  - **GRNBoost2** (Python-based gradient boosting)
  - **ZILGM** (zero-inflated Gaussian graphical model)
  - **PCzinb** (partial correlation with zero-inflated negative binomial model)
  - **JRF** (joint random forests for multi-dataset inference)

### ✨ Major Features

- `infer_networks()`: Infer regulatory networks from count matrices.
- `cutoff_adjacency()`: Apply null model thresholding to weighted adjacency matrices.
- `create_consensus()`: Build consensus networks across methods or datasets.
- `plotROC()`, `community_similarity()`: Performance evaluation tools using ROC curves, AUC, and community structure metrics.
- `plotg()`, `community_path()`: Network visualization functions based on `ggraph`.

### 🧪 Testing and Documentation

- Added runnable examples to all major exported functions.
- Built a comprehensive README and External Dependencies installation guide.
- Set up internal unit tests for core network inference and evaluation functionalities.
- Prepared detailed documentation with reproducible examples.

### 🏛 Project Funding

- Supported by the **National Centre for HPC, Big Data and Quantum Computing** under the European Union – Next Generation EU – CN00000013 (CUP: B93C22000620006).

