#' A wrapper to estimate latent unknown clusters with multi-omics data
#'
#' @param G Required. Genetic features/environmental exposures, a \code{\link{matrix}}.
#' @param Z Required. Biomarkers/other omics data, a \code{\link{matrix}}.
#' @param Y Required. Outcome, it is suggested to transform it into a n by 1 \code{\link{matrix}}.
#' @param CoG Optional, matrix. Covariates to be adjusted for estimating the latent cluster.
#' @param CoY Optional, matrix. Covariates to be adjusted for estimating the outcome.
#' @param family Required. Type of outcome Y. It should be choose from "normal", "binary".
#' @param useY Whether or not to include the information of Y to estimate the latent clusters. Default is TRUE.
#' @param K Required, an array of integers, representing the number of latent clusters
#' @param Rho_G Optional, an array of penalties for variable selection in exposures
#' @param Rho_Z_InvCov Optional, an array of penalties for variable selection in Biomarkers (penalizing inverse of covariance matrix)
#' @param Rho_Z_CovMu Optional, an array of penalties for variable selection in Biomarkers (penalizing mean)
#' @param CV Number of folds in K-fold cross validation
#' @param tol Tolerance of convergence for EM algorithm
#' @return A list:
#' \item{best_model}{the best model among different combination of tuning parameters}
#' \item{tune_list}{a data frame contains combination of tuning parameters and corresponding BIC}
#' \item{res_model}{a list stores LUCID model corresponding to each combination of tuning parameters}
#' @export
#'
lucid <- function(G, 
                  Z, 
                  Y, 
                  CoG = NULL, 
                  CoY = NULL, 
                  family = "normal", 
                  useY = TRUE,
                  K = 2:5,
                  Rho_G = 0, 
                  Rho_Z_Mu = 0,
                  Rho_Z_Cov = 0, 
                  ...){
  # combinations of tuning parameters
  tune_list <- expand.grid(K, Rho_G, Rho_Z_Mu, Rho_Z_Cov)
  colnames(tune_list) <- c("K", "Rho_G", "Rho_Z_Mu", "Rho_Z_Cov")
  m <- nrow(tune_list)
  tune_list$BIC <- rep(0, m)
  # fit models for each combination
  res_model <- vector(mode = "list",
                      length = m)
  for(i in 1:m) {
    fit <- try(est.lucid(G = G, 
                         Z = Z, 
                         Y = Y,
                         CoG = CoG, 
                         CoY = CoY,
                         family = family, 
                         useY = useY,
                         K = tune_list[i, 1],
                         Rho_G = tune_list[i, 2],
                         Rho_Z_Mu = tune_list[i, 3],
                         Rho_Z_Cov = tune_list[i, 4],
                         ...))
    if("try-error" %in% class(fit)) {
      tune_list[i, 5] <- NA
    } else {
      tune_list[i, 5] <- summary(fit$BIC)
    }
    res_model[[i]] <- fit
  }
  x <- min(tune_list[, 5], na.rm = TRUE)
  if(is.na(x)) {
    stop("LUCID model fails to converge given current tuning parameters")
  }
  best_model <- res_model[[which(tune_list[, 5]== x)]]
  return(list(best_model = best_model,
              tune_list = tune_list,
              res_model = res_model))
}