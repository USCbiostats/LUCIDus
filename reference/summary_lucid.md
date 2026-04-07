# Summarize results of LUCID model

Summarize results of LUCID model

## Usage

``` r
summary_lucid(object, boot.se = NULL)
```

## Arguments

- object:

  A LUCID model fitted by [`est_lucid`](est_lucid.md)

- boot.se:

  An object returned by [`boot_lucid`](boot_lucid.md), which contains
  the bootstrap confidence intervals

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

# summarize lucid model
summary_lucid(fit1)

# summarize lucid model with bootstrap CIs
summary_lucid(fit1, boot.se = boot1)
} # }
```
