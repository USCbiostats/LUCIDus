# Changelog

## LUCIDus 2.2.1

CRAN release: 2022-11-08

1.  Update vignette; create package website; update citations.

## LUCIDus 2.2.0

### New features

1.  New function `lucid`: an estimation function integrating model and
    variable selection. We recommend user to use this new function for
    estimating LUCID models. It calls the workhorse functions
    `est_lucid` and `tune_lucid` in the back-end. \## Major changes
2.  Rename `est.lucid` and `boot.lucid` to `est_lucid` and `boot_lucid`.
3.  Update vignette accordingly.
4.  Use `mclust` to choose the optimal geometric model for omics data Z
    by default.

## LUCIDus 2.1.5-2

CRAN release: 2022-06-01

### Minor bugs fixed

1.  Update package dependencies
2.  Fix bugs when running examples

## LUCIDus 2.1.5-1

CRAN release: 2022-04-06

### Minor bugs fixed

1.  Fix bugs when running tests

## LUCIDus 2.1.5

CRAN release: 2022-03-28

### New features

1.  Add progress bar for `boot.lucid` function to track how many
    iterations are done.
2.  Add `verbose` parameter in `est.lucid` to disable automatic output
    in R console.

### Other changes

1.  Fix bug for `boot.lucid`
2.  Fix bug for `predict_lucid`

## LUCIDus 2.1.4

CRAN release: 2022-03-02

### New features

1.  `est.lucid` is enhanced to:

- Deal with missing values in omics data, inclduing sporadict missing
  pattern and list-wise missing pattern.
- Integrated imputation for missing values.
- Different approaches for initialization of EM algorithm, including
  mclust and guess from uniform distribution.

2.  New function `lucid`, a wrapper function to perform model selection
    over grid of K and penalty terms.
3.  `plot.lucid`: new option to change color for nodes and lines
4.  `boot.lucid`: return original output from
    [`boot::boot`](https://rdrr.io/pkg/boot/man/boot.html) to allow user
    plot diagnostic figures.

### Other changes

1.  Use log-sum-exponential trick and update the likelihood function to
    avoid under/overflow
2.  Major bug fixes for `boot.lucid` function
3.  Update simulated dataset
4.  Update corresponding examples in documentation.

## LUCIDus 2.1.0

CRAN release: 2020-07-22

### New features

A new variable selection framework is applied to LUCID. \* A lasso type
penalty is applied on the mean of biomarkers \* A glasso method is
applied on the variance-covariance structure to achieve sparsity
covariance matrix \* We apply a new variable selection criteria, which
takes both mean and coviarnce matrix of biomarkers into account.

### Other changes

- Fix bugs in `pred.lucid()`. Now it can predict both latent cluster and
  the outcome.

## LUCIDus 2.0.0

CRAN release: 2020-05-18

This is a feature update to the whole package. It rewrite all the codes
to make the model fitting procedure much faster (10 to 50 times) than
v1.0.0. Also, the grammar of LUCID changed to a more user-friendly
version. (Please note, this version is not backward compatible)

### New features

- [`est.lucid()`](../reference/est.lucid.md): previously called
  `est_lucid`. Fit the LUCID mode much faster; use mclust to initialize
  and produce a more stable estimate of the model; fix the bugs dealing
  with missing values in biomarker data.
- `summary.lucid()`: previously called `summary_lucid`. An S3 method
  function which can directly be called by `summary`; provide with a
  nice table with detailed interpretation of the model; add option to
  calculate 95% CI based on result returned by
  [`boot.lucid()`](../reference/boot.lucid.md).
- `plot.lucid()`: previously called `plot_lucid`. An S3 method function
  which can be directly called by `plot`; change the color palette.
- `predict.lucid()`: previously called `pred_lucid`. An S3 method
  function which can be directly called by `predict`.
- [`boot.lucid()`](../reference/boot.lucid.md): previously called
  `boot_lucid`. Provide with a neat output.

### Other changes

- Update the vignette.
- Updated the citation after getting published by the
  *[Bioinformatics](https://doi.org/10.1093/bioinformatics/btz667)*;
- Minor bug fixes.
