# LUCIDus 2.0.0


This is a feature update to the whole package. It rewrite all the codes to make the model fitting procedure much faster (10 to 50 times) than v1.0.0. Also, the grammer of LUCID changed to a more user-friendly version. (Please note, this version is not backward compatible)

## New features

* `est.lucid()`: previously called `est_lucid`. Fit the LUCID mode much faster; use mclust to initialize and produce a more stable estimate of the model; fix the bugs dealing with missing values in biomarker data.
* `summary.lucid()`: previously called `summary_lucid`. An S3 method function which can directly be called by `summary`; provide with a nice table with detailed interpretation of the model; add option to calcualte 95% CI based on result returned by `boot.lucid()`.
* `plot.lucid()`: previously called `plot_lucid`. An S3 method function which can be directly called by `plot`; change the color palette.
* `predict.lucid()`: previously called `pred_lucid`. An S3 method function which can be directly called by `predict`.
* `boot.lucid()`: previously called `boot_lucid`. Provide with a neat output.

## Other changes

* Update the vignette.
* Updated the citation after getting published by the *[Bioinformatics](https://doi.org/10.1093/bioinformatics/btz667)*;
* Minor bug fixes.
