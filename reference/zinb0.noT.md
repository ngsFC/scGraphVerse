# Structure learning with zero-inflated negative binomial model (mean only)

This function estimates the adjacency matrix of a ZINB model given a
matrix of counts, using the optim function. Uses BiocParallel for
parallelization.

## Usage

``` r
zinb0.noT(X, maxcard, alpha, extend, nCores = 1)
```

## Arguments

- X:

  the matrix of counts (n times p).

- maxcard:

  the uper bound of the cardinality of the conditional sets K

- alpha:

  the significant level of the tests

- extend:

  if TRUE it considers the union of the tests, otherwise it considers
  the intersection.

- nCores:

  number of cores for parallelization

## Value

the estimated adjacency matrix of the graph.

## Details

This approach assumes that the structure of the graph only depends on
the mean parameter, treating zero inflation as a technical noise effect.
We call this model `zinb0`.
