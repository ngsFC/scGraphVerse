# Infer Gene Regulatory Networks from Expression Matrices

Infers weighted gene regulatory networks (GRNs) from one or more
expression matrices using different inference methods: `"GENIE3"`,
`"GRNBoost2"`, `"ZILGM"`, `"JRF"`, or `"PCzinb"`.

## Usage

``` r
infer_networks(
  count_matrices_list,
  method = c("GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb"),
  adjm = NULL,
  nCores = 1,
  grnboost_modules = NULL,
  genie3_params = list(),
  grnboost2_params = list(),
  zilgm_params = list(),
  jrf_params = list(),
  pczinb_params = list(),
  verbose = FALSE
)
```

## Arguments

- count_matrices_list:

  A MultiAssayExperiment object containing expression data from multiple
  experiments or conditions.

- method:

  Character string. Inference method to use. One of: `"GENIE3"`,
  `"GRNBoost2"`, `"ZILGM"`, `"JRF"`, or `"PCzinb"`.

- adjm:

  Optional. Reference adjacency matrix for matching dimensions when
  using `"ZILGM"` or `"PCzinb"`.

- nCores:

  Integer. Number of CPU cores to use for parallelization. Default: 1.

- grnboost_modules:

  Python modules required for `GRNBoost2` (created via reticulate).

- genie3_params:

  List of parameters for GENIE3 method:

  - `regulators`: Vector of regulator gene names (default: all)

  - `targets`: Vector of target gene names (default: all genes)

  - `treeMethod`: "RF" or "ET" (default: "RF")

  - `K`: Number of candidate regulators (default: "sqrt")

  - `nTrees`: Number of trees per ensemble (default: 1000)

  - `seed`: Random seed for reproducibility (default: NULL)

- grnboost2_params:

  List of parameters for GRNBoost2 method:

  - `tf_names`: Vector of transcription factor names (default:all)

  - `gene_names`: Vector of target gene names (default: all)

  - `client_or_address`: Dask client or address (default: NULL)

  - `seed`: Random seed for reproducibility (default: NULL)

- zilgm_params:

  List of parameters for ZILGM method:

  - `lambda`: Regularization parameter (default: 0.1)

  - `alpha`: Elastic net mixing parameter (default: 1)

  - `max_iter`: Maximum iterations (default: 100)

  - `tol`: Convergence tolerance (default: 1e-4)

- jrf_params:

  List of parameters for JRF method:

  - `ntree`: Number of trees (default: 1000)

  - `mtry`: Number of variables to sample at each split (default:
    sqrt(p))

- pczinb_params:

  List of parameters for PCzinb method:

  - `gamma`: Regularization parameter (default: 0.1)

  - `beta`: Beta parameter (default: 0.1)

  - `max_iter`: Maximum iterations (default: 100)

  - `tol`: Convergence tolerance (default: 1e-4)

- verbose:

  Logical. If TRUE, display messages. Default: FALSE.

## Value

A list of inferred networks:

- For `"GENIE3"`, `"GRNBoost2"`, `"ZILGM"`, and `"PCzinb"`, a list of
  inferred network objects (edge lists or adjacency matrices).

- For `"JRF"`, a list of data frames with inferred edge lists for each
  condition or dataset.

## Details

Each expression matrix is preprocessed automatically depending on its
object type (`Seurat`, `SingleCellExperiment`, or plain matrix).

Parallelization behavior:

- **GENIE3**: No external parallelization; internal `nCores` parameter
  controls computation.

- **ZILGM**: Uses `nCores` parameter for internal parallelization.

- **GRNBoost2** and **PCzinb**: Parallelized across matrices using
  BiocParallel.

- **JRF**: Joint modeling of all matrices together using optimized C
  implementation.

Methods are based on:

- **GENIE3**: Random Forest-based inference (Huynh-Thu et al., 2010).

- **GRNBoost2**: Gradient boosting trees using arboreto (Moerman et al.,
  2019).

- **ZILGM**: Zero-Inflated Graphical Models for scRNA-seq (Zhang et al.,
  2021).

- **JRF**: Joint Random Forests across multiple conditions (Petralia et
  al., 2015).

- **PCzinb**: Pairwise correlation under ZINB models (Nguyen et al.,
  2023).

## Examples

``` r
data("toy_counts")

# Infer networks (toy_counts is already a MultiAssayExperiment)
networks <- infer_networks(
    count_matrices_list = toy_counts,
    method = "GENIE3",
    nCores = 1
)
head(networks[[1]])
#>   regulatoryGene targetGene    weight
#> 1          HLA-B        FTL 0.1959007
#> 2           CD74      CXCR4 0.1554003
#> 3          HLA-A      HLA-B 0.1522556
#> 4          HLA-B      HLA-A 0.1519417
#> 5            FTL       FTH1 0.1411901
#> 6           FTH1        FTL 0.1288961
```
