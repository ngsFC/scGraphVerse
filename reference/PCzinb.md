# Structure learning for count data using PC algorithms

This function performs structure learning for count data using various
PC algorithms adapted for different distributional assumptions including
Poisson, Negative Binomial, and Zero-Inflated Negative Binomial models.

## Usage

``` r
PCzinb(
  x,
  method = c("poi", "nb", "zinb0", "zinb1"),
  alpha = NULL,
  maxcard = 2,
  extend = TRUE,
  nCores = 1,
  whichAssay = "processed",
  ...
)
```

## Arguments

- x:

  A matrix of count data (n × p), SummarizedExperiment, or
  SingleCellExperiment object. For matrix input, rows are samples and
  columns are genes.

- method:

  The algorithm used to estimate the graph: `poi` (Poisson), `nb`
  (Negative Binomial), `zinb0` (Zero-Inflated NB with structure only in
  mu), or `zinb1` (Zero-Inflated NB with structure in both mu and pi).

- alpha:

  The significance level of the tests. Default: 2 \* pnorm(nrow(x)^0.2,
  lower.tail = FALSE).

- maxcard:

  The upper bound of the cardinality of the conditional sets K. Default:
  2.

- extend:

  If TRUE, considers the union of the tests; if FALSE, considers the
  intersection. Default: TRUE.

- nCores:

  Number of cores for parallel processing. Default: 1.

- whichAssay:

  The assay to use as input (for SummarizedExperiment or
  SingleCellExperiment objects). Default: "processed".

- ...:

  Additional arguments (currently unused).

## Value

- If x is a matrix: the estimated adjacency matrix of the graph

- If x is a SummarizedExperiment: the object with adjacency matrix
  stored in metadata as `adj_mat`

- If x is a SingleCellExperiment: the object with adjacency matrix
  stored as rowPair

## Details

PCzinb performs structure learning using PC algorithms for count data.
Different methods handle different distributional assumptions:

- `poi`: Poisson distribution

- `nb`: Negative Binomial distribution

- `zinb0`: Zero-Inflated NB with structure only in mean parameter

- `zinb1`: Zero-Inflated NB with structure in both mean and
  zero-inflation parameters

For SummarizedExperiment and SingleCellExperiment inputs, if the
specified `whichAssay` is "processed" but not found, the function will
use the first assay and issue a warning recommending QPtransform().

## Examples

``` r
# Matrix input
mat <- matrix(rpois(50, 5), nrow = 10)
PCzinb(mat, method = "poi")
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    0    0    0    0
#> [2,]    0    0    0    1    0
#> [3,]    0    0    0    0    0
#> [4,]    0    1    0    0    0
#> [5,]    0    0    0    0    0

# SummarizedExperiment input
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
se <- SummarizedExperiment(matrix(rpois(50, 5), ncol = 10))
se_result <- PCzinb(se, method = "poi")
#> Warning: We recommend to use QPtransform() before learning the graph.

# SingleCellExperiment input
library(SingleCellExperiment)
sce <- SingleCellExperiment(matrix(rpois(50, 5), ncol = 10))
sce_result <- PCzinb(sce, method = "poi")
#> Warning: We recommend to use QPtransform() before learning the graph.
rowPair(sce_result)
#> SelfHits object with 2 hits and 1 metadata column:
#>            from        to |         x
#>       <integer> <integer> | <numeric>
#>   [1]         3         5 |         1
#>   [2]         5         3 |         1
#>   -------
#>   nnode: 5
```
