# (NB) model. The NB distribution is parametrized by two parameters: the mean value and the dispersion of the negative binomial distribution

(NB) model. The NB distribution is parametrized by two parameters: the
mean value and the dispersion of the negative binomial distribution

## Usage

``` r
nb.OptimizeDispersion(mu, Y, n)
```

## Arguments

- mu:

  the vector mean of the negative binomial

- Y:

  the vector of counts

- n:

  the length of the vector to return Note that theta is sometimes called
  inverse dispersion parameter (and phi=1/theta is then called the
  dispersion parameter). We follow the convention that the variance of
  the NB variable with mean mu and dispersion theta is mu + mu^2/theta.

## Value

A vector of length n with the optimized dispersion parameter values.
