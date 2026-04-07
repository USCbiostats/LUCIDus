# Impute missing data by optimizing the likelihood function

Impute missing data by optimizing the likelihood function

## Usage

``` r
fill_data(obs, mu, sigma, p, index)
```

## Arguments

- obs:

  a vector of length M

- mu:

  a matrix of size M x K

- sigma:

  a matrix of size M x M x K

- p:

  a vector of length K

- index:

  a vector of length M, indicating whether a value is missing or not in
  the raw data

## Value

an observation with updated imputed value
