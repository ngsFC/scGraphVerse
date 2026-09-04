# Compute Network Topological Properties

Calculates topological properties (modularity, number of communities,
density, transitivity) for one or more networks.

## Usage

``` r
compute_topology_metrics(control_output, predicted_list)
```

## Arguments

- control_output:

  A list output from
  [`community_path()`](https://ngsfc.github.io/scGraphVerse/reference/community_path.md)
  representing the ground truth network. Must contain a `graph` (igraph
  object) and `communities$membership`.

- predicted_list:

  A list of lists, each output from
  [`community_path()`](https://ngsfc.github.io/scGraphVerse/reference/community_path.md)
  representing predicted networks.

## Value

A list containing:

- `topology_measures`: A data frame with Modularity,Communities,
  Density, and Transitivity for each prediction.

- `control_topology`: A numeric vector of the same metrics for the
  control network.

## Details

This function requires the **igraph** package. Topological metrics
include:

- Modularity: Quality of community division

- Communities: Number of detected communities

- Density: Proportion of actual vs. possible edges

- Transitivity: Clustering coefficient

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
comm_truth <- community_path(toy_adj_matrix)
#> Detecting communities...

#> Running pathway enrichment...
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> Warning: qvalue::qvalue() failed, returning NA for qvalue. Error: missing values and NaN's not allowed if 'na.rm' is FALSE
#> 'select()' returned 1:1 mapping between keys and columns

topo <- compute_topology_metrics(comm_truth, list(comm_cons))
```
