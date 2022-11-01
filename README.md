
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LUCIDus: Integreted clustering with multi-view data

<!-- badges: start -->

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/LUCIDus?color=green)](https://cran.r-project.org/package=LUCIDus)
![](https://cranlogs.r-pkg.org/badges/grand-total/LUCIDus?color=blue)
[![](https://raw.githubusercontent.com/USCbiostats/badges/master/tommy-image-badge.svg)](https://image.usc.edu)
<!-- badges: end -->

The **LUCIDus** package implements the statistical method LUCID proposed
in the research paper [A Latent Unknown Clustering Integrating
Multi-Omics Data (LUCID) with Phenotypic
Traits](https://doi.org/10.1093/bioinformatics/btz667)
(*Bioinformatics*, 2020). LUCID conducts integrated clustering by using
multi-view data, including exposures, omics data with/without outcome.
**LUCIDus** features variable selection, incorporating missingness in
omics data, visualization of LUCID model via Sankey diagram, bootstrap
inference and functions for tuning model parameters.

If you are interested in integration of omic data to estiamte mediator
or latent structures, please check out [Conti
Lab](https://contilab.usc.edu/about/) to learn more.

## Installation

You can install the development version of LUCIDus from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("USCbiostats/LUCIDus")

# install package with vignettes
devtools::install_github("USCbiostats/LUCIDus", build_vignettes = TRUE)
```

## Usage

Please refer to the
[tutorial](https://yinqi93.github.io/LUCIDus/articles/LUCID_vignette.html).

## Acknowledgments

-   Dr. Lida Chatzi and awesome researchers from [Chatzi
    Lab](https://twitter.com/chatzilab?lang=en)
-   USC IMAGE Group (Supported by the National Cancer Institute at the
    National Institutes of Health Grant P01 CA196569)
-   Dr. Zhao Yang
