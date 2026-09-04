# Classify Edges as TP, FP, or FN

Compares a consensus network to a reference network and classifies edges
as True Positives (TP), False Positives (FP), or False Negatives (FN).

## Usage

``` r
classify_edges(consensus_matrix, reference_matrix = NULL, use_stringdb = TRUE)
```

## Arguments

- consensus_matrix:

  A
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object representing the consensus network.

- reference_matrix:

  A
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object representing the reference (ground truth) network. If `NULL`,
  STRINGdb is used to generate a high-confidence physical network.

- use_stringdb:

  Logical. If TRUE and `reference_matrix` is NULL, queries STRINGdb for
  reference network. Default: TRUE.

## Value

A list containing:

- `tp_edges`: Character vector of True Positive edges

- `fp_edges`: Character vector of False Positive edges

- `fn_edges`: Character vector of False Negative edges

- `consensus_graph`: igraph object of consensus network

- `reference_graph`: igraph object of reference network

- `edge_colors`: Color vector for TP (red) and FN (blue) edges

- `use_stringdb`: Logical indicating if STRINGdb was used

## Details

If `reference_matrix` is NULL and `use_stringdb = TRUE`, this function
queries STRINGdb to generate a human high-confidence (score \> 900)
physical interaction network.

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

# classify_edges expects SummarizedExperiment objects
edge_class <- classify_edges(consensus, ref_se)
```
