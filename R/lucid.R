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
#'
#' @return A list:
#' \item best_model: the best model among different combination of tuning parameters
#' \item tune_list: a data frame contains combination of tuning parameters and corresponding BIC
#' \item res_model: a list stores LUCID model corresponding to each combination of tuning parameters
#' @export
#'
#' @examples
lucid <- function(G, Z, Y, CoG = NULL, CoY = NULL, family = "normal", useY = TRUE,
                  K = 2:5, Rho_G = 0, Rho_Z_InvCov = 0, Rho_Z_CovMu = 0){
  # combinations of tuning parameters
  tune_list <- expand.grid(K, Rho_G, Rho_Z_CovMu, Rho_Z_InvCov)
  colnames(tune_list) <- c("K", "Rho_G", "Rho_Z_InvCov", "Rho_Z_CovMu")
  m <- nrow(tune_list)
  tune_list$BIC <- rep(NA, m)
  # fit models for each combination
  res_model <- lapply(1:m, function(x){
    K <- tune_list[x, 1]
    Rho_G <- tune_list[x, 2]
    Rho_Z_InvCov <- tune_list[x, 3]
    Rho_Z_CovMu <- tune_list[x, 4]
    selectG <- ifelse(Rho_G != 0, TRUE, FALSE)
    selectZ <- ifelse(Rho_Z_CovMu == 0 && Rho_Z_InvCov == 0, FALSE, TRUE)
    fit <- try(est.lucid(G = G, Z = Z, Y = Y, CoG = CoG, CoY = CoY,
                         family = family, useY = useY, K = K,
                         tune = def.tune(Rho_G = Rho_G,
                                         Rho_Z_CovMu = Rho_Z_CovMu,
                                         Rho_Z_InvCov = Rho_Z_InvCov,
                                         Select_Z = selectZ,
                                         Select_G = selectG)))
    return(fit)
  })
  # choose the best model
  for (i in 1:m) {
    fit <- res_model[[i]]
    if(!("try-error" %in% class(fit))){
      tune_list$BIC[i] <- summary(fit)$BIC
    }
  }
  x <- min(tune_list$BIC, na.rm = TRUE)
  best_model <- res_model[[which(tune_list$BIC == x)]]
  return(list(best_model = best_model,
              tune_list = tune_list,
              res_model = res_model))
}