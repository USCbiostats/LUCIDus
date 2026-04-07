# A simulated dataset for LUCID

This is an example dataset to illustrate LUCID model. It is simulated by
assuming there are 2 latent clusters in the data. We assume the
exposures are associated with latent cluster which ultimately affects
the PFAS concentration and liver injury in children. The latent clusters
are also characterized by differential levels of metabolites.

## Usage

``` r
sim_data
```

## Format

A list with 5 matrices corresponding to exposures (G), omics data (Z), a
continuous outcome, a binary outcome and 2 covariates (can be used
either as CoX or CoY). Each matrice contains 2000 observations.

- G:

  10 exposures

- Z:

  10 metabolites

- Y_normal:

  Outcome, PFAS concentration in children

- Y_bninary:

  Bianry outcome, liver injury status

- Covariates:

  2 continous covariates, can be treated as either CoX or CoY

- X:

  Latent clusters
