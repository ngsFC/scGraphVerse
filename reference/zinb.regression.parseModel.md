# Parse ZINB regression model

Given the parameters of a ZINB regression model, this function parses
the model and computes the vector of log(mu), logit(pi), and the
dimensions of the different components of the vector of parameters.

## Usage

``` r
zinb.regression.parseModel(alpha, A.mu, A.pi)
```

## Arguments

- alpha:

  the vectors of parameters c(a.mu, a.pi) concatenated

- A.mu:

  matrix of the model (default=empty)

- A.pi:

  matrix of the model (default=empty)

## Value

A list with slots `logMu`, `logitPi`, `dim.alpha` (a vector of length 2
with the dimension of each of the vectors `a.mu`, `a.pi` in `alpha`),
and `start.alpha` (a vector of length 2 with the starting indices of the
2 vectors in `alpha`)
