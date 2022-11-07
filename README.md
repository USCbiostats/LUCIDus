
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

You can install the stable release from CRAN by:

``` r
install.packages("LUCIDus")
```

or install the development version of LUCIDus from
[GitHub](https://github.com/USCbiostats/LUCIDus) with:

``` r
# install.packages("devtools")
devtools::install_github("USCbiostats/LUCIDus")

# install package with vignettes
devtools::install_github("USCbiostats/LUCIDus", build_vignettes = TRUE)
```

## Usage

Please refer to the
[tutorial](https://USCbiostats.github.io/LUCIDus/articles/LUCIDus.html).

## Citation

    #> 
    #> To cite LUCID methods, please use:
    #> 
    #>   Cheng Peng, Jun Wang, Isaac Asante, Stan Louie, Ran Jin, Lida Chatzi,
    #>   Graham Casey, Duncan C Thomas, David V Conti (2019). A latent unknown
    #>   clustering integrating multi-omics data (LUCID) with phenotypic
    #>   traits. Bioinformatics, btz667. URL
    #>   https://doi.org/10.1093/bioinformatics/btz667
    #> 
    #> To cite LUCIDus R package, please use:
    #> 
    #>   Yinqi Zhao (2022). LUCIDus: an R package to implement the LUCID
    #>   model. CRAN. R package version 2.2 URL
    #>   https://yinqi93.github.io/LUCIDus/index.html
    #> 
    #> To see these entries in BibTeX format, use 'print(<citation>,
    #> bibtex=TRUE)', 'toBibtex(.)', or set
    #> 'options(citation.bibtex.max=999)'.

## Acknowledgments

-   Dr. Lida Chatzi and awesome researchers from [Chatzi
    Lab](https://twitter.com/chatzilab?lang=en)
-   [Dr. Nikos
    Stratakis](https://www.isglobal.org/en/our-team/-/profiles/29105)
-   Dr. Cheng Peng
-   [USC IMAGE Group](https://p01.uscbiostatistics.org/) (Supported by
    the National Cancer Institute at the National Institutes of Health
    Grant P01 CA196569)
