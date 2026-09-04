# Compute Community Assignment Similarity Metrics

Calculates community assignment similarity between a reference community
structure and one or more predicted structures using variation of
information (VI), normalized mutual information (NMI), and adjusted Rand
index (ARI).

## Usage

``` r
compute_community_metrics(control_output, predicted_list)
```

## Arguments

- control_output:

  A list output from
  [`community_path()`](https://ngsfc.github.io/scGraphVerse/reference/community_path.md)
  representing the ground truth network. Must contain
  `communities$membership`.

- predicted_list:

  A list of lists, each output from
  [`community_path()`](https://ngsfc.github.io/scGraphVerse/reference/community_path.md)
  representing predicted networks to compare.

## Value

A data frame with columns VI, NMI, and ARI for each prediction. Row
names indicate which prediction (e.g., "Predicted_1").

## Details

This function requires the **igraph** package. Lower VI values indicate
better similarity (VI = 0 is perfect match). Higher NMI and ARI values
indicate better similarity (both range 0-1).

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
    nCores = 1
)

consensus <- create_consensus(binary_se, method = "union")
comm_cons <- community_path(consensus)
#> Detecting communities...

#> Running pathway enrichment...
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
comm_truth <- community_path(toy_adj_matrix)
#> Detecting communities...

#> Running pathway enrichment...
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> Warning: qvalue::qvalue() failed, returning NA for qvalue. Error: missing values and NaN's not allowed if 'na.rm' is FALSE
#> 'select()' returned 1:1 mapping between keys and columns

metrics <- compute_community_metrics(comm_truth, list(comm_cons))
```
