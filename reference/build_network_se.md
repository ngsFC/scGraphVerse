# Create a SummarizedExperiment for Network Storage

Constructs a SummarizedExperiment container for multiple gene regulatory
network adjacency matrices with shared gene space (p×p matrices).

## Usage

``` r
build_network_se(
  networks,
  networkData = NULL,
  geneData = NULL,
  metadata = list()
)
```

## Arguments

- networks:

  A list of adjacency matrices (all must have same dimensions)

- networkData:

  A DataFrame with metadata for each network

- geneData:

  Optional. A DataFrame with gene-level annotations

- metadata:

  Optional. List of global metadata

## Value

A SummarizedExperiment object where each assay is a p×p network

## Examples

``` r
# Example 1: Building SE from a list of adjacency matrices
data("toy_adj_matrix")

# Create a list of network matrices
net_list <- list(
    network1 = toy_adj_matrix,
    network2 = toy_adj_matrix
)

# Build SummarizedExperiment
network_se <- build_network_se(net_list)
network_se
#> class: SummarizedExperiment 
#> dim: 35 35 
#> metadata(1): object_type
#> assays(2): network1 network2
#> rownames(35): ACTG1 ARPC2 ... UBA52 UBC
#> rowData names(1): gene
#> colnames(35): ACTG1 ARPC2 ... UBA52 UBC
#> colData names(1): gene

# Example 2: Using with inferred networks
data("toy_counts")

networks <- infer_networks(
    count_matrices_list = toy_counts,
    method = "GENIE3",
    nCores = 1
)

# generate_adjacency() internally uses build_network_se()
wadj_se <- generate_adjacency(networks)
wadj_se
#> class: SummarizedExperiment 
#> dim: 35 35 
#> metadata(2): type object_type
#> assays(3): network_1 network_2 network_3
#> rownames(35): ACTG1 ARPC2 ... UBA52 UBC
#> rowData names(1): gene
#> colnames(35): ACTG1 ARPC2 ... UBA52 UBC
#> colData names(1): gene
```
