# Parse ZINB regression model

Given the parameters of a NB regression model, this function parses the
model and computes the vector of log(mu), and the dimensions of the
different components of the vector of parameters. See
[`nb.loglik.regression`](https://ngsfc.github.io/scGraphVerse/reference/nb.loglik.regression.md)
for details of the NB regression model and its parameters.

## Usage

``` r
nb.regression.parseModel(alpha, A.mu)
```

## Arguments

- alpha:

  the vectors of parameters c(a.mu) concatenated

- A.mu:

  matrix of the model (default=empty)

## Value

A list with slot `logMu`,

## See also

[`nb.loglik.regression`](https://ngsfc.github.io/scGraphVerse/reference/nb.loglik.regression.md)
