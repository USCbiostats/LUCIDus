% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_lucid.R
\name{summary_lucid}
\alias{summary_lucid}
\title{Summarize results of LUCID model}
\usage{
summary_lucid(object, boot.se = NULL)
}
\arguments{
\item{object}{A LUCID model fitted by \code{\link{est_lucid}}}

\item{boot.se}{An object returned by \code{\link{boot_lucid}}, 
which contains the bootstrap confidence intervals}
}
\description{
Summarize results of LUCID model
}
\examples{
\dontrun{
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
}
}
