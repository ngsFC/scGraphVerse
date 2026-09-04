# Gradient of the log-likelihood of the NB regression model

This function computes the gradient of the log-likelihood of a NB
regression model given a vector of counts.

## Usage

``` r
nb.loglik.regression.gradient(
  alpha,
  Y,
  A.mu = matrix(nrow = length(Y), ncol = 0),
  C.theta = matrix(0, nrow = length(Y), ncol = 1)
)
```

## Arguments

- alpha:

  the vectors of parameters a.mu concatenated

- Y:

  the vector of counts

- A.mu:

  matrix of the model (see Details, default=empty)

- C.theta:

  matrix of the model (see Details, default=zero)

## Value

The gradient of the log-likelihood.

## Details

The regression model is described in
[`nb.loglik.regression`](https://ngsfc.github.io/scGraphVerse/reference/nb.loglik.regression.md).

## See also

[`nb.loglik.regression`](https://ngsfc.github.io/scGraphVerse/reference/nb.loglik.regression.md)
