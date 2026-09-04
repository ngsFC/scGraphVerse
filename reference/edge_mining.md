# Edge Mining of Gene Interactions Using PubMed

Query PubMed for literature evidence supporting predicted gene–gene
interactions.

## Usage

``` r
edge_mining(
  predicted_list,
  ground_truth,
  delay = 1,
  query_field = "Title/Abstract",
  query_edge_types = c("TP", "FP", "FN"),
  max_retries = 10,
  BPPARAM = BiocParallel::bpparam()
)
```

## Arguments

- predicted_list:

  A list of predicted adjacency matrices (row and column names are gene
  symbols), or a
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object containing adjacency matrices.

- ground_truth:

  A 0/1 adjacency matrix with row and column names.

- delay:

  Numeric. Seconds to wait between consecutive queries (default = 1).

- query_field:

  Character. PubMed search field. Options: "Title/Abstract" (default),
  "Title", "Abstract".

- query_edge_types:

  Character vector. Edge types to query: c("TP", "FP", "FN") (default
  all).

- max_retries:

  Integer. Max retries for PubMed queries (default = 10).

- BPPARAM:

  A BiocParallel parameter object. Default = bpparam().

## Value

A named list of data.frames. Each data.frame has columns:

- gene1:

  First gene in interaction

- gene2:

  Second gene

- edge_type:

  One of "TP", "FP", or "FN"

- pubmed_hits:

  Number of PubMed hits

- PMIDs:

  Comma-separated PubMed IDs or NA

- query_status:

  One of "hits_found", "no_hits", or "error"

## Details

This function compares predicted adjacency matrices against a ground
truth matrix, identifies edge types (TP, FP, FN), and queries PubMed for
each gene pair. Returns counts of hits, PMIDs, and query status.

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
#> 1          HLA-B        FTL 0.1960129
#> 2          HLA-A      HLA-B 0.1652298
#> 3           FTH1        FTL 0.1436177
#> 4           CD74      CXCR4 0.1434554
#> 5          HLA-B      HLA-A 0.1432814
#> 6            FTL       FTH1 0.1372178

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
#> [Method: GENIE3] Matrix 1 → Cutoff = 0.06688
#> [Method: GENIE3] Matrix 2 → Cutoff = 0.06206
#> [Method: GENIE3] Matrix 3 → Cutoff = 0.06642
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
em <- edge_mining(consensus, toy_adj_matrix, query_edge_types = "TP")
```
