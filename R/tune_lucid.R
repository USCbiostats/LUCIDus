#' Grid search for tuning parameters to fit the LUCID model
#'
#' @param G Genetic features/environmental exposures, a \code{\link{matrix}}.
#' @param Z Biomarkers/other omics data, a \code{\link{matrix}}.
#' @param Y Disease outcome, it is suggested to transform it into a n by 1 \code{\link{matrix}}.
#' @param CoG Optional, matrix. Covariates to be adjusted for estimating the latent cluster.
#' @param CoY Optional, matrix. Covariates to be adjusted for estimating the outcome.
#' @param K Numeric sequence. Number of latent clusters.
#' @param family Type of outcome Y. It should be choose from "normal", "binary".
#' @param useY Whether or not to include the information of Y to estimate the latent clusters. Default is TRUE.
#' @param Rho_G Numeric sequence, Lasso type penalty for selection of G.
#' @param Rho_Z_InvCov Numeric sequence, Lasso type penalty for the inverse covariance structure of Z.
#' @param Rho_Z_CovMu Numeric sequence, Lasso type penalty for the product of covariance matrix and mean of Z
#'
#' @return A list. Containing model BICs of different combination of tuning parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' tuenpar <- tune.lucid(G = G1, Z = Z1, Y = Y1, family = "binary",
#' Rho_G = seq(0.01, 0.02, by = 0.005),
#' Rho_Z_InvCov = seq(0.1, 0.3, by = 0.1),
#' Rho_Z_CovMu = seq(80, 100, by = 10))
#' }
tune.lucid <- function(G, Z, Y, CoG = NULL, CoY = NULL, family = "normal", useY = TRUE,
                       K = 2:6, Rho_G = 0, Rho_Z_InvCov = 0, Rho_Z_CovMu = 0){
  res <- data.frame(K = K)
  bic <- NULL
  if(length(K) != 1){
    for (i in 1:length(K)) {
      invisible(capture.output(temp.fit <- est.lucid(G = G, Z = Z, Y = Y, CoG = CoG, CoY = CoY,
                                                     family = family, useY = useY, K = K[i])))
      bic <- c(bic, summary(temp.fit)$BIC)
    }
    res$BIC <- bic
    opt.K <- K[bic = which.min(bic)]
    opt.tune <- c(rep(NA, 3), opt.K, min(bic))
    names(opt.tune) <- c("Rho_G", "Rho_Z_InvCov", "Rho_Z_CovMu", "K", "BIC")
  } else{
    opt.K = K
  }

  res2 <- NULL
  if(!is.null(Rho_G) && !is.null(Rho_Z_CovMu) && !is.null(Rho_Z_InvCov)){
    for (i in 1:length(Rho_G)) {
      for(j in 1:length(Rho_Z_InvCov)){
        for (k in 1:length(Rho_Z_CovMu)) {
          invisible(capture.output(temp.fit <- est.lucid(G = G, Z = Z, Y = Y, CoG = CoG, CoY = CoY,
                                                         family = family, useY = useY, K = opt.K,
                                                         tune = def.tune(Rho_G = Rho_G[i], Select_G = TRUE,
                                                                         Rho_Z_CovMu = Rho_Z_CovMu[k],
                                                                         Rho_Z_InvCov = Rho_Z_InvCov[j],
                                                                         Select_Z = TRUE))))
          aa <- c(Rho_G[i], Rho_Z_InvCov[j], Rho_Z_CovMu[k], opt.K, summary(temp.fit)$BIC)
          res2 <- rbind(res2, aa)
        }
      }
    }
    res2 <- as.data.frame(res2)
    colnames(res2) <- c("Rho_G", "Rho_Z_InvCov", "Rho_Z_CovMu", "K", "BIC")
    row.names(res2) <- 1:nrow(res2)
    opt.tune <- res2[which.min(res2$BIC), ]
  }
  return(list(res.K = res,
              res.tune = res2,
              optimal = opt.tune))
}