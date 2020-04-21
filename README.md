
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Introduction

The LUCIDus R package is an integrative tool to obtain a joint
estimation of latent or unknown clusters/subgroups with multi-omics data
and phenotypic traits. This package is an implementation for the novel
statistical method proposed in the research paper “[A Latent Unknown
Clustering Integrating Multi-Omics Data (LUCID) with Phenotypic
Traits](https://doi.org/10.1093/bioinformatics/btz667)” published by the
*Bioinformatics*.

## Citation

Cheng Peng, Jun Wang, Isaac Asante, Stan Louie, Ran Jin, Lida Chatzi,
Graham Casey, Duncan C Thomas, David V Conti, A Latent Unknown
Clustering Integrating Multi-Omics Data (LUCID) with Phenotypic Traits,
Bioinformatics, , btz667,
<https://doi.org/10.1093/bioinformatics/btz667>

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# devtools::install_github("Yinqi93/LUCIDus")
```

## Example

``` r
library(LUCIDus2)
```

The two main functions: `est.lucid()` and `boot.lucid()` are used for
model fitting and estimation of SE of the model parameters. You can also
achieve variable selection by setting tuning parameters in `def.lucid`.
The model outputs can be summarized and visualized using `summary` and
`plot` respectively. Predictions could be made with `pred`.

Estimating latent clusters with multi-omics data, missing values in
biomarker data are allowed, and information in the outcome of interest
can be integrated. For illustration, we use a testing dataset with 10
genetic features (5 causal) and 10 biomarkers (5 causal)

### Integrative clustering without feature selection

First, fit the model with `est.lucid`.

``` r
set.seed(10)
myfit <- est.lucid(G = G1,Z = Z1,Y = Y1, CoY = CovY, K = 2, family = "binary")
myfit
```

![](man/figures/fit1.png) Check the model features.

``` r
summary(myfit)
```

A summary of results start with this: ![](man/figures/sum1.png) Then
visualize the results with Sankey diagram using `plot_lucid()`

``` r
plot(myfit)
```

![](man/figures/Sankey1.png)

### Integrative clustering with feature selection

Run LUCID with tuning parameters and select informative features

``` r
set.seed(10)
myfit2 <- est.lucid(G = G1, Z = Z1, Y = Y1, CoY = CovY, K = 2, family = "binary", useY = FALSE, tune = def.tune(Select_Z = TRUE, Rho_Z_InvCov = 0.2, Rho_Z_CovMu = 90, Select_G = TRUE, Rho_G = 0.02))
selectG <- myfit2$select$selectG
selectZ <- myfit2$select$selectZ
```

Re-fit with selected features

``` r
set.seed(10)
myfit3 <- est.lucid(G = G1[, selectG], Z = Z1[, selectZ], Y = Y1, CoY = CovY, K = 2, family = "binary", useY = FALSE)
```

``` r
plot(myfit3)
```

![](man/figures/Sankey2.png)

### Bootstrap method to obtain SEs for LUCID parameter estimates

``` r
set.seed(10)
myboot <- boot.lucid(G = G1[, selectG], Z = Z1[, selectZ], Y = Y1, CoY = CovY, model = myfit3, R = 50)
summary(myfit3, boot.se = myboot)
```

A detailed summary with 95% CI is provided as below.
![](man/figures/sum2.png)

For more details, see documentations for each function in the R package.

## Built With

  - [devtools](https://cran.r-project.org/package=devtools) - Tools to
    Make Developing R Packages Easier
  - [roxygen2](https://cran.r-project.org/package=roxygen2) - In-Line
    Documentation for R

## Versioning

The current version is 1.0.0.

For the versions available, see the
[Release](https://github.com/Yinqi93/LUCIDus/releases) on this
repository.

## Authors

  - **Yinqi Zhao**

## License

This project is licensed under the GPL-2 License.

## Acknowledgments

  - Cheng Peng, Ph.D.
  - David V. Conti, Ph.D.
  - Zhao Yang, Ph.D.
  - USC IMAGE P1 Group
