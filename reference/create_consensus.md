# Create a Consensus Adjacency Matrix from Multiple Networks

Builds a consensus adjacency matrix from networks stored in a
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
using one of three methods: `"vote"`, `"union"`, or `"INet"`.

## Usage

``` r
create_consensus(
  adj_matrix_list,
  method = c("vote", "union", "INet"),
  weighted_list = NULL,
  theta = 0.04,
  threshold = 0.5,
  ncores = 1,
  tolerance = 0.1,
  nitermax = 50,
  verbose = FALSE
)
```

## Arguments

- adj_matrix_list:

  A
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object containing binary adjacency matrices (square, 0/1) with
  identical dimensions and matching row/column names, or a list of such
  matrices.

- method:

  Character string specifying the consensus strategy. One of:

  - `"vote"` (default): An edge is included if supported by at least
    `threshold` fraction of matrices.

  - `"union"`: An edge is included if present in any matrix.

  - `"INet"`: Combines normalized weighted matrices using
    [`consensusNet`](https://rdrr.io/pkg/INetTool/man/consensusNet.html).

- weighted_list:

  A
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object containing weighted adjacency matrices (required if
  `method = "INet"`), or a list of such matrices.

- theta:

  Numeric. Tuning parameter passed to `consensusNet` (default: `0.04`).

- threshold:

  Numeric between 0 and 1. Threshold for "vote" and "INet" methods.
  Default is `0.5`.

- ncores:

  Integer. Number of CPU cores to use when `method = "INet"`. Default is
  `1`.

- tolerance:

  Numeric. Tolerance for differences between similar graphs in INet
  method. Default is `0.1`.

- nitermax:

  Integer. Maximum number of iterations for INet algorithm. Default is
  `50`.

- verbose:

  Logical. If TRUE, display verbose output for INet method. Default is
  `FALSE`.

## Value

A
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with a single assay containing the consensus adjacency matrix
(binary or weighted, depending on the method). Metadata includes
consensus method and parameters.

## Details

Consensus construction depends on the selected method:

- **vote**:

  Counts the presence of each edge across all matrices and includes
  edges supported by at least `threshold × N` matrices.

- **union**:

  Includes any edge that appears in any matrix.

- **INet**:

  Multiplies binary matrices by corresponding weighted matrices,
  normalizes the results, and applies `consensusNet` to generate a
  consensus network.

For "INet", both binary and weighted adjacency matrices must be provided
with matching dimensions.

## Examples

``` r
data(toy_counts)


# Infer networks (toy_counts is already a MultiAssayExperiment)
networks <- infer_networks(
    count_matrices_list = toy_counts,
    method = "GENIE3",
    nCores = 1
)
head(networks[[1]])
#>   regulatoryGene targetGene    weight
#> 1          HLA-B        FTL 0.2118853
#> 2          HLA-A      HLA-B 0.1642582
#> 3            FTL       FTH1 0.1544068
#> 4          HLA-B      HLA-A 0.1424878
#> 5           CD74      CXCR4 0.1415349
#> 6           FTH1        FTL 0.1392659

# Generate adjacency matrices
wadj_se <- generate_adjacency(networks)
swadj_se <- symmetrize(wadj_se, weight_function = "mean")

# Apply cutoff
binary_se <- cutoff_adjacency(
    count_matrices = toy_counts,
    weighted_adjm_list = swadj_se,
    n = 1,
    method = "GENIE3",
    quantile_threshold = 0.95,
    nCores = 1,
    debug = TRUE
)
#> [Method: GENIE3] Matrix 1 → Cutoff = 0.06341
#> [Method: GENIE3] Matrix 2 → Cutoff = 0.06468
#> [Method: GENIE3] Matrix 3 → Cutoff = 0.06116
head(binary_se[[1]])
#> [1] "ACTG1" "ARPC2" "ARPC3" "BTF3"  "CD3D"  "CD3E" 

consensus <- create_consensus(binary_se, method = "union")
head(consensus)
#> class: SummarizedExperiment 
#> dim: 6 35 
#> metadata(4): type method threshold object_type
#> assays(1): consensus
#> rownames(6): ACTG1 ARPC2 ... CD3D CD3E
#> rowData names(1): gene
#> colnames(35): ACTG1 ARPC2 ... UBA52 UBC
#> colData names(1): gene
```
