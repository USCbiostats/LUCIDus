# Inference of LUCID model based on bootstrap resampling

Generate `R` bootstrap replicates of LUCID parameters and derive
confidence interval (CI) base on bootstrap. Bootstrap replicates are
generated based on nonparameteric resampling, implemented by `ordinary`
method of codeboot::boot function.

## Usage

``` r
boot_lucid(G, Z, Y, CoG = NULL, CoY = NULL, model, conf = 0.95, R = 100)
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

## Value

A list, containing the following components:

- beta:

  effect estimate for each exposure

- mu:

  cluster-specific mean for each omics feature

- gamma:

  effect estiamte for the association btween latent cluster and outcome

- bootstrap:

  The `boot` object returned by `boot:boot`

## Examples

``` r
if (FALSE) { # \dontrun{
# use simulated data
G <- sim_data$G
Z <- sim_data$Z
Y_normal <- sim_data$Y_normal

# fit lucid model
fit1 <- est_lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, 
seed = 1008)

# conduct bootstrap resampling
boot1 <- boot_lucid(G = G, Z = Z, Y = Y_normal, model = fit1, R = 100)

# check distribution for bootstrap replicates of the variable of interest
plot(boot1$bootstrap, 1)

# use 90% CI
boot2 <- boot_lucid(G = G, Z = Z, Y = Y_normal, model = fit1, R = 100, conf = 0.9)
} # }
```
