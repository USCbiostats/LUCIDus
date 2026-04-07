# Deprecated function boot.lucid

This function deprecates. Please use boot_lucid instead.

## Usage

``` r
boot.lucid(G, Z, Y, CoG = NULL, CoY = NULL, model, conf = 0.95, R = 100)
```

## Arguments

- G:

  Exposures, a numeric vector, matrix, or data frame. Categorical
  variable should be transformed into dummy variables. If a matrix or
  data frame, rows represent observations and columns correspond to
  variables.

- Z:

  Omics data, a numeric matrix or data frame. Rows correspond to
  observations and columns correspond to variables.

- Y:

  Outcome, a numeric vector. Categorical variable is not allowed. Binary
  outcome should be coded as 0 and 1.

- CoG:

  Optional, covariates to be adjusted for estimating the latent cluster.
  A numeric vector, matrix or data frame. Categorical variable should be
  transformed into dummy variables.

- CoY:

  Optional, covariates to be adjusted for estimating the association
  between latent cluster and the outcome. A numeric vector, matrix or
  data frame. Categorical variable should be transformed into dummy
  variables.

- model:

  A LUCID model fitted by `est.lucid`.

- conf:

  A numeric scalar between 0 and 1 to specify confidence level(s) of the
  required interval(s).

- R:

  An integer to specify number of bootstrap replicates for LUCID model.
  If feasible, it is recommended to set R \>= 1000.
