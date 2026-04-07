# I-step of LUCID

Impute missing data in Z by maximizing the likelihood given fixed
parameters of LUCID

## Usage

``` r
Istep_Z(Z, p, mu, sigma, index)
```

## Arguments

- Z:

  an N by P matrix representing the omics data

- p:

  an N by K matrix representing posterior inclusion probability for each
  latent cluster

- mu:

  an M by K matrix representing cluster-specific means

- sigma:

  an M by M by K array representing cluster-specific covariance

- index:

  an N by M matrix representing missing values in Z

## Value

a complete dataset of Z
