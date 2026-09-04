# Generate Adjacency Matrices from Gene Interaction Tables

Constructs adjacency matrices from a list of data frames (network edge
lists) and returns them in a
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object.

## Usage

``` r
generate_adjacency(df_list, nCores = 1)
```

## Arguments

- df_list:

  A list of data frames. Each data frame must have three columns:

  Gene1

  :   Character. First gene in the interaction.

  Gene2

  :   Character. Second gene in the interaction.

  Weight

  :   Numeric. Weight or strength of the interaction from `Gene1` to
      `Gene2`.

- nCores:

  Integer. Number of CPU cores to use for parallel processing. Defaults
  to the number of available workers from the current BiocParallel
  backend.

## Value

A
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object where each assay is a square numeric adjacency matrix (p×p
genes). Diagonal entries are set to zero (no self-interactions).

## Details

The function first identifies all unique genes across all data frames to
define the matrix dimensions. For each interaction table, it places the
corresponding weights at the appropriate gene-pair positions.
Parallelization is handled by BiocParallel for improved performance on
multiple datasets.

Missing weights (`NA`) are ignored during construction. Only gene pairs
matching the global gene list are inserted.

## Examples

``` r
data("toy_counts")

# Infer networks (toy_counts is already a MultiAssayExperiment)
networks <- infer_networks(
    count_matrices_list = toy_counts,
    method = "GENIE3",
    nCores = 1
)
head(networks[[1]])
#>   regulatoryGene targetGene    weight
#> 1          HLA-B        FTL 0.1941378
#> 2           CD74      CXCR4 0.1579812
#> 3          HLA-A      HLA-B 0.1559264
#> 4            FTL       FTH1 0.1533053
#> 5          HLA-B      HLA-A 0.1470368
#> 6           FTH1        FTL 0.1400307

# Generate adjacency matrices
wadj_se <- generate_adjacency(networks) # returns SummarizedExperiment
head(wadj_se[[1]])
#> [1] "ACTG1" "ARPC2" "ARPC3" "BTF3"  "CD3D"  "CD3E" 
```
