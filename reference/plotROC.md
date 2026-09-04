# Plot ROC Curves for Inferred Networks

Computes and visualizes Receiver Operating Characteristic (ROC) curves
for predicted adjacency matrices stored in a
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object, compared against a binary ground truth network.

## Usage

``` r
plotROC(predicted_se, ground_truth, plot_title, is_binary = FALSE)
```

## Arguments

- predicted_se:

  A
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object containing predicted adjacency matrices as assays. Each matrix
  must share dimnames with `ground_truth`;entries may be binary (0/1) or
  continuous weights.

- ground_truth:

  A square binary matrix indicating true interactions (1) in the upper
  triangle. Must match dims and names of each assay in `predicted_se`.

- plot_title:

  Character string. Title for the ROC plot.

- is_binary:

  Logical. If `TRUE`, treat matrices as binary predictions. Default
  `FALSE` for weighted predictions.

## Value

A list with:  
`auc`: data frame of AUC per matrix.  
`plot`: the ROC plot (via ggplot2).

## Details

For binary matrices, a single TPR/FPR point is computed per matrix. For
weighted ones, a full ROC curve is computed from continuous scores.
Diagonals are ignored; symmetry is not enforced.

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
#> 1          HLA-B        FTL 0.1983718
#> 2          HLA-A      HLA-B 0.1613441
#> 3            FTL       FTH1 0.1465988
#> 4           CD74      CXCR4 0.1427142
#> 5          HLA-B      HLA-A 0.1389132
#> 6           FTH1        FTL 0.1340791

# Generate and symmetrize adjacency matrices (returns SummarizedExperiment)
wadj_se <- generate_adjacency(networks)
swadj_se <- symmetrize(wadj_se, weight_function = "mean")

# plotROC now accepts SummarizedExperiment directly
roc_res <- plotROC(
    swadj_se,
    toy_adj_matrix,
    plot_title = "ROC Curve: GENIE3",
    is_binary = FALSE
)
#> Registered S3 methods overwritten by 'pROC':
#>   method    from
#>   print.roc fmsb
#>   plot.roc  fmsb
roc_res$plot

auc_joint <- roc_res$auc
```
