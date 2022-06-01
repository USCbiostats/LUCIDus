# LUCIDus 2.1.5-2
## Minor bugs fixed
1. Update package dependencies
2. Fix bugs when running examples

# LUCIDus 2.1.5-1
## Minor bugs fixed
1. Fix bugs when running tests

# LUCIDus 2.1.5
## New features
1. Add progress bar for `boot.lucid` function to track how many iterations are done.
2. Add `verbose` parameter in `est.lucid` to disable automatic output in R console.

## Other changes
1. Fix bug for `boot.lucid`
2. Fix bug for `predict_lucid`

# LUCIDus 2.1.4
## New features
1. `est.lucid` is enhanced to:
  * Deal with missing values in omics data, inclduing sporadict missing pattern and list-wise missing pattern.
  * Integrated imputation for missing values.
  * Different approaches for initialization of EM algorithm, including mclust and guess from uniform distribution.
2. New function `lucid`, a wrapper function to perform model selection over grid of K and penalty terms.
3. `plot.lucid`: new option to change color for nodes and lines
4. `boot.lucid`: return original output from `boot::boot` to allow user plot diagnostic figures.

## Other changes
1. Use log-sum-exponential trick and update the likelihood function to avoid under/overflow
2. Major bug fixes for `boot.lucid` function
3. Update simulated dataset 
3. Update corresponding examples in documentation.


# LUCIDus 2.1.0
## New features
A new variable selection framework is applied to LUCID. 
* A lasso type penalty is applied on the mean of biomarkers
* A glasso method is applied on the variance-covariance structure to achieve sparsity covariance matrix
* We apply a new variable selection criteria, which takes both mean and coviarnce matrix of biomarkers into account.

## Other changes
* Fix bugs in `pred.lucid()`. Now it can predict both latent cluster and the outcome.



# LUCIDus 2.0.0


This is a feature update to the whole package. It rewrite all the codes to make the model fitting procedure much faster (10 to 50 times) than v1.0.0. Also, the grammar of LUCID changed to a more user-friendly version. (Please note, this version is not backward compatible)

## New features

* `est.lucid()`: previously called `est_lucid`. Fit the LUCID mode much faster; use mclust to initialize and produce a more stable estimate of the model; fix the bugs dealing with missing values in biomarker data.
* `summary.lucid()`: previously called `summary_lucid`. An S3 method function which can directly be called by `summary`; provide with a nice table with detailed interpretation of the model; add option to calculate 95% CI based on result returned by `boot.lucid()`.
* `plot.lucid()`: previously called `plot_lucid`. An S3 method function which can be directly called by `plot`; change the color palette.
* `predict.lucid()`: previously called `pred_lucid`. An S3 method function which can be directly called by `predict`.
* `boot.lucid()`: previously called `boot_lucid`. Provide with a neat output.

## Other changes

* Update the vignette.
* Updated the citation after getting published by the *[Bioinformatics](https://doi.org/10.1093/bioinformatics/btz667)*;
* Minor bug fixes.
