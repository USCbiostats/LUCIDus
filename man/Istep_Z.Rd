% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missing.R
\name{Istep_Z}
\alias{Istep_Z}
\title{I-step of LUCID}
\usage{
Istep_Z(Z, p, mu, sigma, index)
}
\arguments{
\item{Z}{an N by P matrix representing the omics data}

\item{p}{an N by K matrix representing posterior inclusion probability for each 
latent cluster}

\item{mu}{an M by K matrix representing cluster-specific means}

\item{sigma}{an M by M by K array representing cluster-specific covariance}

\item{index}{an N by M matrix representing missing values in Z}
}
\value{
a complete dataset of Z
}
\description{
Impute missing data in Z by maximizing the likelihood given fixed
parameters of LUCID
}
