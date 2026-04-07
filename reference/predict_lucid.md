# Predict cluster assignment and outcome based on LUCID model

Predict cluster assignment and outcome based on LUCID model

## Usage

``` r
predict_lucid(model, G, Z, Y = NULL, CoG = NULL, CoY = NULL, response = TRUE)
```

## Arguments

- model:

  A model fitted and returned by [`est_lucid`](est_lucid.md)

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

- response:

  If TRUE, when predicting binary outcome, the response will be
  returned. If FALSE, the linear predictor is returned.

## Value

A list contains predicted latent cluster and outcome for each
observation

## Examples

``` r
if (FALSE) { # \dontrun{
# prepare data
G <- sim_data$G
Z <- sim_data$Z
Y_normal <- sim_data$Y_normal

# fit lucid model
fit1 <- est_lucid(G = G, Z = Z, Y = Y_normal, K = 2, family = "normal")

# prediction on training set
pred1 <- predict_lucid(model = fit1, G = G, Z = Z, Y = Y_normal)
pred2 <- predict_lucid(model = fit1, G = G, Z = Z)
} # }
```
