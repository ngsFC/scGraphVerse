# Log-likelihood of negative binomial model, for a fixed dispersion parameter

Given a unique dispersion parameter and a set of counts, together with a
corresponding set of mean parameters, this function computes the sum of
the log-probabilities of the counts under the NB model. The dispersion
parameter is provided to the function through zeta = log(theta), where
theta is sometimes called the inverse dispersion parameter.

## Usage

``` r
nb.loglik.dispersion(zeta, Y, mu)
```

## Arguments

- zeta:

  a vector, the log of the inverse dispersion parameters of the negative
  binomial model

- Y:

  a vector of counts

- mu:

  a vector of mean parameters of the negative binomial

## Value

the log-likelihood of the model.
