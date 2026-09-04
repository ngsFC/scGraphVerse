# Visualize Graphs from Adjacency Matrices

Generates and arranges multiple graph visualizations from a list of
adjacency matrices. Each matrix is converted into an undirected igraph
object and visualized using a force-directed layout via ggraph.

## Usage

``` r
plotg(adj_matrix_list)
```

## Arguments

- adj_matrix_list:

  A list of square, symmetric adjacency matrices with zeros on the
  diagonal (no self-loops), or a
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object containing such matrices as assays. Each matrix represents an
  undirected graph.

## Value

A grid of plots displaying all valid graphs in the input list.

## Details

Each adjacency matrix is validated to ensure it is square and symmetric.
Disconnected nodes (degree zero) are removed prior to visualization.
Graphs are visualized with a force-directed layout using ggraph, and
multiple plots are arranged into a grid with gridExtra.

Each subplot title includes the graph index, number of nodes, and number
of edges.

## Note

This function requires the following packages: igraph, ggraph, and
gridExtra. If any are missing, an informative error will be thrown.

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
#> 1          HLA-B        FTL 0.2047071
#> 2          HLA-A      HLA-B 0.1535271
#> 3           CD74      CXCR4 0.1515372
#> 4          HLA-B      HLA-A 0.1405707
#> 5           FTH1        FTL 0.1345138
#> 6            FTL       FTH1 0.1307083

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
#> [Method: GENIE3] Matrix 1 → Cutoff = 0.06565
#> [Method: GENIE3] Matrix 2 → Cutoff = 0.06893
#> [Method: GENIE3] Matrix 3 → Cutoff = 0.06881
head(binary_se[[1]])
#> [1] "ACTG1" "ARPC2" "ARPC3" "BTF3"  "CD3D"  "CD3E" 
plotg(binary_se)
```
