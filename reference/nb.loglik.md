# Log-likelihood of the negative binomial model Given a vector of counts, this function computes the sum of the log-probabilities of the counts under a negative binomial (NB) model. The NB distribution is parametrized by two parameters: the mean value and the dispersion of the negative binomial distribution

Log-likelihood of the negative binomial model Given a vector of counts,
this function computes the sum of the log-probabilities of the counts
under a negative binomial (NB) model. The NB distribution is
parametrized by two parameters: the mean value and the dispersion of the
negative binomial distribution

## Usage

``` r
nb.loglik(Y, mu, theta)
```

## Arguments

- Y:

  the vector of counts

- mu:

  the vector of mean parameters of the negative binomial

- theta:

  the vector of dispersion parameters of the negative binomial, or a
  single scalar is also possible if the dispersion parameter is
  constant. Note that theta is sometimes called inverse dispersion
  parameter (and phi=1/theta is then called the dispersion parameter).
  We follow the convention that the variance of the NB variable with
  mean mu and dispersion theta is mu + mu^2/theta.

## Value

the log-likelihood of the model.
