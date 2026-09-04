# Compare Consensus and Reference Graphs or STRINGdb Networks

Convenience wrapper that classifies edges and visualizes the comparison
between consensus and reference networks. For more control, use the
individual functions:
[`classify_edges`](https://ngsfc.github.io/scGraphVerse/reference/classify_edges.md)
and
[`plot_network_comparison`](https://ngsfc.github.io/scGraphVerse/reference/plot_network_comparison.md).

## Usage

``` r
compare_consensus(
  consensus_matrix,
  reference_matrix = NULL,
  false_plot = FALSE
)
```

## Arguments

- consensus_matrix:

  A
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object representing the consensus network.

- reference_matrix:

  Optional. A
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  obj representing the reference network. If `NULL`, STRINGdb is
  queried.

- false_plot:

  Logical. If `TRUE`, displays False Positives plot. Default is `FALSE`.

## Value

A `ggplot` object visualizing the comparison.

## Details

If no `reference_matrix` is provided, STRINGdb is queried to generate a
high-confidence physical interaction network.

## Note

Requires ggraph and ggplot2. If `reference_matrix` is NULL, also
requires STRINGdb. If `false_plot = TRUE`, requires patchwork.

## See also

[`classify_edges`](https://ngsfc.github.io/scGraphVerse/reference/classify_edges.md),
[`plot_network_comparison`](https://ngsfc.github.io/scGraphVerse/reference/plot_network_comparison.md)

## Examples

``` r
data(toy_counts)
data(toy_adj_matrix)


# Infer networks (toy_counts is already a MultiAssayExperiment)
networks <- infer_networks(
    count_matrices_list = toy_counts,
    method = "GENIE3",
    nCores = 1
)
head(networks[[1]])
#>   regulatoryGene targetGene    weight
#> 1          HLA-B        FTL 0.2160796
#> 2          HLA-A      HLA-B 0.1684571
#> 3           CD74      CXCR4 0.1654010
#> 4          HLA-B      HLA-A 0.1529067
#> 5           FTH1        FTL 0.1462899
#> 6            FTL       FTH1 0.1377266

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
#> [Method: GENIE3] Matrix 1 → Cutoff = 0.06609
#> [Method: GENIE3] Matrix 2 → Cutoff = 0.06367
#> [Method: GENIE3] Matrix 3 → Cutoff = 0.07024
head(binary_se[[1]])
#> [1] "ACTG1" "ARPC2" "ARPC3" "BTF3"  "CD3D"  "CD3E" 

consensus <- create_consensus(binary_se, method = "union")

# Wrap reference matrix in SummarizedExperiment
ref_se <- build_network_se(list(reference = toy_adj_matrix))

# Compare consensus to reference
compare_consensus(
    consensus,
    reference_matrix = ref_se,
    false_plot = FALSE
)

```
