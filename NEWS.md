# scGraphVerse News

## scGraphVerse 0.99.21 (2025-10-21)

### Bug Fixes
- Added `doRNG` to Imports to support GENIE3 parallel processing with `nCores > 1`
- Fixed igraph deprecation warning by using `mode = "max"` instead of `mode = "undirected"` for asymmetric matrices
- Fixed ggplot2 deprecation warning by replacing `size` with `linewidth` in `geom_line()` calls
- Added `BiocStyle` to Suggests for vignette building
- Added `useDynLib(scGraphVerse, .registration = TRUE)` to NAMESPACE for proper C code registration

### Improvements
- Enhanced NAMESPACE imports to satisfy R CMD check requirements
- Improved vignette dependency management

## scGraphVerse 0.99.0 (Development Version)

### 🚀 Enhanced Parameter Support & Internal Dependencies

#### New Features
- **Enhanced Parameter Support**: All main functions now expose the complete parameter sets from underlying tools:
  - `infer_networks()`: Now supports all method-specific parameters (GENIE3, GRNBoost2, ZILGM, JRF, PCzinb)
  - `create_consensus()`: Now supports all INetTool parameters (tolerance, nitermax, verbose)
  - `community_path()`: Now supports all robin package parameters (method-specific args, comparison params, parallel processing)
- **Verbose Output**: Added verbose parameter to track progress across all functions
- **Seed Support**: Added seed parameter for reproducible results
- **Internal Dependencies**: Integrated ZILGM and JRF functions internally with proper GPL-2 attribution

#### Internal Dependencies Integration
- **ZILGM Functions**: Integrated `find_lammax()` and `zilgm()` functions from ZILGM package (Park et al., 2021)
- **JRF Functions**: Added fallback implementation for JRF functionality when package not available
- **INet-Tool**: Maintained as external dependency (already in imports)
- **License Update**: Changed to GPL (>= 2) to accommodate integrated GPL-2 code

#### Breaking Changes
- License changed from MIT to GPL (>= 2) due to integrated GPL-2 code
- ZILGM no longer required as external dependency (functions integrated)
- JRF now optional with internal fallback implementation

#### Bug Fixes & Improvements
- Improved parameter validation and error handling
- Better documentation with comprehensive parameter descriptions
- Enhanced examples showing new parameter capabilities

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

