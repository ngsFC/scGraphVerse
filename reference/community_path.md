# Community Detection and Pathway Enrichment Analysis

Detects gene communities within an adjacency network using one or two
community detection methods, and performs pathway enrichment for each
detected community.

## Usage

``` r
community_path(
  adj_matrix,
  methods = "louvain",
  pathway_db = "KEGG",
  organism = c("human", "mouse"),
  genes_path = 5,
  plot = TRUE,
  verbose = TRUE,
  method_params = list(),
  comparison_params = list(),
  BPPARAM = BiocParallel::bpparam()
)
```

## Arguments

- adj_matrix:

  A square adjacency matrix or a
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object containing a single adjacency matrix. Row and column names must
  correspond to gene symbols.

- methods:

  A character vector of one or two community detection methods supported
  by robin. If two are given, performance is compared and the best is
  selected. Default: `"louvain"`.

- pathway_db:

  Character string specifying the pathway database to use: `"KEGG"` or
  `"Reactome"`. Default: `"KEGG"`.

- organism:

  Character string specifying the organism: `"human"` or `"mouse"`.
  Default: `"human"`.

- genes_path:

  Integer. Minimum number of genes per community to run enrichment
  analysis. Default: `5`.

- plot:

  Logical. If `TRUE`, generates a plot of detected communities. Default:
  `TRUE`.

- verbose:

  Logical. If `TRUE`, shows progress messages. Default: `TRUE`.

- method_params:

  List of parameters for community detection methods. Common parameters
  include:

  - `resolution`: Resolution parameter for Louvain/Leiden (default: 1)

  - `steps`: Number of steps for Walktrap (default: 4)

  - `spins`: Number of spins for Spinglass (default: 25)

  - `nb.trials`: Number of trials for Infomap (default: 10)

- comparison_params:

  List of parameters for robin comparison:

  - `measure`: Stability measure ("vi", "nmi", "split.join",
    "adjusted.rand"). Default: "vi"

  - `type`: Robin construction type ("dependent", "independent").
    Default: "independent"

  - `rewire.w.type`: Rewiring strategy for weighted graphs. Default:
    "Rewire"

- BPPARAM:

  BiocParallel backend for parallel processing. Default:
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).

## Value

A list with elements:

- `communities`: List with `best_method` and a named vector of community
  membership per gene.

- `pathways`: List of enrichment results per community (only for
  communities meeting size threshold).

- `graph`: The igraph object with community annotations.

## Details

If two methods are provided, the function uses `robinCompare` and
selects the method with higher AUC. Pathway enrichment is done via
clusterProfiler (KEGG) or via ReactomePA (Reactome). Communities smaller
than `genes_path` are excluded.

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
#> 1          HLA-B        FTL 0.2045461
#> 2           CD74      CXCR4 0.1579789
#> 3          HLA-A      HLA-B 0.1547400
#> 4          HLA-B      HLA-A 0.1519506
#> 5            FTL       FTH1 0.1423699
#> 6           FTH1        FTL 0.1396498

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
    nCores = 1,
    debug = TRUE
)
#> [Method: GENIE3] Matrix 1 → Cutoff = 0.06602
#> [Method: GENIE3] Matrix 2 → Cutoff = 0.06785
#> [Method: GENIE3] Matrix 3 → Cutoff = 0.06985

# Create consensus (returns SummarizedExperiment)
consensus <- create_consensus(binary_se, method = "union")

# community_path now accepts SummarizedExperiment objects directly
comm_cons <- community_path(consensus)
#> 
#> 
#> Detecting communities...

#> Running pathway enrichment...
#> 'select()' returned 1:1 mapping between keys and columns
#> Reading KEGG annotation online: "https://rest.kegg.jp/link/hsa/pathway"...
#> Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/hsa"...
#> kegg_category.rda is not found, download it online...
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
```
