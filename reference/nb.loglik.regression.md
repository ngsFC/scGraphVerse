# log-likelihood of the NB regression model

This function computes the log-likelihood of a NB regression model given
a vector of counts.

## Usage

``` r
nb.loglik.regression(
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

  matrix of the model (\\log(\theta)\\, default=zero)

## Value

the log-likelihood.

## Details

The regression model is parametrized as follows: \$\$log(\mu) = A\_\mu
\* a\_\mu\$\$ \$\$log(\theta) = C\_\theta\$\$ where \\\mu, \theta\\ are
respectively the vector of mean parameters of the NB distribution, and
the vector of inverse dispersion parameters. The log-likelihood of a
vector of parameters \\\alpha = a\_\mu\\
