# Simulate Zero-Inflated Negative Binomial (ZINB) Count Matrices with Sequencing Depth

Simulates one or more count matrices following a zero-inflated negative
binomial (ZINB) distribution, incorporating gene-gene interaction
structures and cell-specific sequencing depth variation.

## Usage

``` r
zinb_simdata(
  n,
  p,
  B,
  mu_range,
  mu_noise,
  theta,
  pi,
  kmat = 1,
  depth_range = NA
)
```

## Arguments

- n:

  Integer. Number of cells (samples) in each simulated matrix.

- p:

  Integer. Number of genes (features) in each simulated matrix.

- B:

  A symmetric binary adjacency matrix (0/1) defining gene-gene
  connectivity. Row and column names correspond to gene names.

- mu_range:

  List of numeric vectors (length 2 each). Range of gene expression
  means for each simulated matrix.

- mu_noise:

  Numeric vector. Mean of background noise for each matrix.

- theta:

  Numeric vector. Dispersion parameters of the negative binomial
  distribution for each matrix. Smaller `theta` implies higher
  overdispersion.

- pi:

  Numeric vector. Probability of excess zeros (`0 < pi < 1`) for each
  matrix.

- kmat:

  Integer. Number of count matrices to simulate. Default is 1.

- depth_range:

  Numeric vector of length 2 or `NA`. Range of total sequencing depth
  per cell. If `NA`, no depth adjustment is performed.

## Value

A list containing `kmat` matrices. Each matrix has:

- Rows representing cells (`cell_1, ..., cell_n`).

- Columns representing genes (`rownames(B)`).

- Count values following a ZINB distribution.

## Details

Each simulated matrix:

1.  Generates gene expression values based on a ZINB model.

2.  Modulates expression using the adjacency matrix `B`.

3.  Applies random sequencing depth scaling if `depth_range` is
    provided.

Useful for benchmarking single-cell RNA-seq network inference methods
with dropout events and network structure.

## Examples

``` r
data(toy_adj_matrix)
nodes <- nrow(toy_adj_matrix)
sims <- zinb_simdata(
    n = 50,
    p = nodes,
    B = toy_adj_matrix,
    mu_range = list(c(1, 4), c(1, 7), c(1, 10)),
    mu_noise = c(1, 3, 5),
    theta = c(1, 0.7, 0.5),
    pi = c(0.2, 0.2, 0.2),
    kmat = 3,
    depth_range = c(0.8 * nodes * 3, 1.2 * nodes * 3)
)
```
