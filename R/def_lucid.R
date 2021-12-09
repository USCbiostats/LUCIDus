#' Control parameters for EM algorithm
#' 
#' @param tol Convergence criteria for the EM algorithm. Default is 0.001.
#' @param max_itr Maximum number of iterations in each try of fitting process, integer, default is 1000.
#' @param max_tot.itr Maximum number of total iterations, integer, default is 10000.
#'
#' @export

def.control <- function(tol = 1e-3,
                        max_itr = 1000,
                        max_tot.itr = 10000) {
  if(tol > 1e-2) {
    stop("The convergence tolerance is too large, please use a smaller value than 1e-2")
  }
  
  list(tol = tol,
       max_itr = max_itr,
       max_tot.itr = max_tot.itr)
}



#' Define tuning parameters of regularization for LUCID model.
#'
#' @param Rho_G Numeric. Penalty for selecting exposures G.
#' @param Rho_Z_InvCov Numeric. Penalty to obtain a sparse covariance matrix for omics data Z.
#' @param Rho_Z_Mu Numeric. 
#' @param Rho_Z_CovMu Deprecated argument. Penalty for the product of the inverse of the covariance of omicss, which will produce a sparse matrix for the mean.
#' @param Select_G Deprecated argument. Flag for variable selection in genetic features/environmental exposures.
#' @param Select_Z Deprecated argument. Flag for variable selection in omics data Z.
#'
#' @export

def.tune <- function(Rho_G = 0,
                     Rho_Z_InvCov = 0,
                     Rho_Z_Mu = 0,
                     Rho_Z_CovMu = 0,
                     Select_G = FALSE,
                     Select_Z = FALSE) {
  if(!missing(Select_G) | !missing(Select_Z)) {
    warning("arguments 'Select_G' and 'Select_Z' are deprecated")
  }
  
  if(Select_G) {
    Rho_G == Rho_G
  }
  if(Select_Z) {
    Rho_Z_InvCov = Rho_Z_InvCov
    Rho_Z_CovMu = Rho_Z_CovMu
  }
  structure(list(Select_G = Select_G,
                 Select_Z = Select_Z,
                 Rho_G = Rho_G,
                 Rho_Z_InvCov = Rho_Z_InvCov,
                 Rho_Z_CovMu = Rho_Z_CovMu))
}



