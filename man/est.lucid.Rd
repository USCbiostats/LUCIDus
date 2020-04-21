% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_lucid.R
\name{est.lucid}
\alias{est.lucid}
\title{Estimate latent unknown clusters with multi-omics data}
\usage{
est.lucid(
  G,
  Z,
  Y,
  CoG = NULL,
  CoY = NULL,
  K = 2,
  family = "normal",
  useY = TRUE,
  control = def.control(),
  tune = def.tune(),
  Z.var.str = NULL
)
}
\arguments{
\item{G}{Genetic features/environmental exposures, a \code{\link{matrix}}.}

\item{Z}{Biomarkers/other omics data, a \code{\link{matrix}}.}

\item{Y}{Disease outcome, it is suggested to transform it into a n by 1 \code{\link{matrix}}.}

\item{CoG}{Optional, matrix. Covariates to be adjusted for estimating the latent cluster.}

\item{CoY}{Optional, matrix. Covariates to be adjusted for estimating the outcome.}

\item{K}{Number of latent clusters.}

\item{family}{Type of outcome Y. It should be choose from "normal", "binary".}

\item{useY}{Whether or not to include the information of Y to estimate the latent clusters. Default is TRUE.}

\item{control}{A list of tolerance parameters used by EM algorithm. See \code{\link{def_control}}.}

\item{tune}{A list of tuning parameters used by variable selection procedure. See \code{\link{def_tune}}}

\item{Z.var.str}{The variance-covariance structure for the biomarkers. See \code{\link{mclustModelNames}} for details.}
}
\value{
A list which contains the several features of LUCID, including:
\item{pars}{Estimates of parameters of LUCID, including beta (estimates of genetic feature/environmental exposure), mu (estimates of cluster-specific biomarker means), sigma (estimates of the cluster-specific biomarker variance-covariance matrix) and gamma(estimates of cluster-specific effect and covariates effect related to the outcome)}
\item{K}{Number of latent cluster}
\item{Z.var.str}{The model used to estimate the cluster-specific variance-covariance matrix, for further details, see \code{\link{mclust}}}
\item{likelihood}{The log likelihood of the LUCID model}
\item{post.p}{Predicted probability of belonging to each latent cluster}
}
\description{
This function estimates the latent clusters by integrating genetic features/environmental exposures, biomarkers with/without the outcome of interest. Variable selection is available for analyzing the high-dimensional data.
}
\examples{
set.seed(10)
fit1 <- est.lucid(G = G1,Z = Z1,Y = Y1, CoY = CovY, K = 2, family = "binary")
fit2 <- est.lucid(G = G1,Z = Z1,Y = Y1, CoY = CovY, K = 2, family = "binary", tune = def.tune(Select_Z = TRUE, Rho_Z_InvCov = 0.1, Rho_Z_CovMu = 90, Select_G = TRUE, Rho_G = 0.02))
}
\references{
Cheng Peng, Jun Wang, Isaac Asante, Stan Louie, Ran Jin, Lida Chatzi, Graham Casey, Duncan C Thomas, David V Conti, A Latent Unknown Clustering Integrating Multi-Omics Data (LUCID) with Phenotypic Traits, Bioinformatics, , btz667, https://doi.org/10.1093/bioinformatics/btz667.
}
\author{
Yinqi Zhao, Cheng Peng, Zhao Yang, David V. Conti
}