# Visualize Network Comparison

Creates visualization plots comparing consensus and reference networks,
showing True Positives (TP), False Negatives (FN), and optionally False
Positives (FP) edges.

## Usage

``` r
plot_network_comparison(edge_classification, show_fp = FALSE)
```

## Arguments

- edge_classification:

  A list output from
  [`classify_edges()`](https://ngsfc.github.io/scGraphVerse/reference/classify_edges.md)
  containing edge classifications and graph objects.

- show_fp:

  Logical. If TRUE, displays False Positive edges in a separate plot.
  Default: FALSE.

## Value

A `ggplot` object visualizing the comparison. If `show_fp = TRUE`, a
combined plot using patchwork is returned.

## Details

This function requires the **ggraph** and **ggplot2** packages. If
`show_fp = TRUE`, the **patchwork** package is also required.

The plots differentiate:

- TP/CE (True Positives/Confirmed Edges): Red edges present in both

- FN/ME (False Negatives/Missing Edges): Blue edges in reference only

- FP/EE (False Positives/Extra Edges): Edges in consensus only

If STRINGdb was used for reference, labels are CE/ME/EE. Otherwise,
TP/FN/FP.

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

# Generate and symmetrize adjacency matrices (returns SummarizedExperiment)
wadj_se <- generate_adjacency(networks)
swadj_se <- symmetrize(wadj_se, weight_function = "mean")

# Apply cutoff (returns SummarizedExperiment)
binary_se <- cutoff_adjacency(
    count_matrices = toy_counts,
    weighted_adjm_list = swadj_se,
    n = 1,
    method = "GENIE3",
    quantile_threshold = 0.95,
    nCores = 1
)

# Create consensus (returns SummarizedExperiment)
consensus <- create_consensus(binary_se, method = "union")

# Wrap reference matrix in SummarizedExperiment
ref_se <- build_network_se(list(reference = toy_adj_matrix))

# Classify edges
edge_class <- classify_edges(consensus, ref_se)

# Plot comparison
plot_network_comparison(edge_class, show_fp = FALSE)
```
