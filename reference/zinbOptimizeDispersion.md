# (ZINB) model. The ZINB distribution is parametrized by three parameters: the mean value and the dispersion of the negative binomial distribution, and the probability of the zero component.

(ZINB) model. The ZINB distribution is parametrized by three parameters:
the mean value and the dispersion of the negative binomial distribution,
and the probability of the zero component.

## Usage

``` r
zinbOptimizeDispersion(mu, logitPi, Y, n)
```

## Arguments

- mu:

  the vector mean of the negative binomial

- logitPi:

  the vector of logit of the probabilities of the zero component

- Y:

  the vector of counts

- n:

  length of the returned vector

## Value

A vector of length n with the optimized dispersion parameter values.
