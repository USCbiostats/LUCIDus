#' Control parameters for EM algorithm
#' 
#' @param tol Convergence criteria for the EM algorithm. Default is 0.001.
#' @param max_itr Maximum number of iterations in each try of fitting process, integer, default is 1000.
#' @param max_tot.itr Maximum number of total iterations, integer, default is 10000.
#'
#' @return A list of tolerance settings for LUCID.
#' @export
#' @author Yinqi Zhao, Cheng Peng, Zhao Yang, David V. Conti

def.control <- function(tol = 1e-3,
                        max_itr = 1000,
                        max_tot.itr = 10000){
  structure(list(tol = tol,
                 max_itr = max_itr,
                 max_tot.itr = max_tot.itr))
}



#' Define tuning parameters of regularization for LUCID model.
#'
#' @param Rho_G Numeric. Penalty for selection on genetic features/environmental exposures.
#' @param Rho_Z_InvCov Numeric. Penalty for the inverse of the covariance of biomarkers, which will produce a sparse matrix.
#' @param Rho_Z_CovMu Numeric. Penalty for the product of the inverse of the covariance of biomarkers, which will produce a sparse matrix for the mean.
#' @param Select_G Flag for variable selection in genetic features/environmental exposures. Default is FALSE.
#' @param Select_Z Flag for variable selection in biomarkers. Default is FALSE.
#'
#' @return A list of tuning parameters and settings will be returned for integrative clustering.
#' @export
#' @author Yinqi Zhao, Cheng Peng, Zhao Yang, David V. Conti

def.tune <- function(Rho_G = 0,
                     Rho_Z_InvCov = 0,
                     Rho_Z_CovMu = 0,
                     Select_G = FALSE,
                     Select_Z = FALSE) {
  if(Select_G){
    Rho_G == Rho_G
  }
  if(Select_Z){
    Rho_Z_InvCov = Rho_Z_InvCov
    Rho_Z_CovMu = Rho_Z_CovMu
  }
  structure(list(Select_G = Select_G,
                 Select_Z = Select_Z,
                 Rho_G = Rho_G,
                 Rho_Z_InvCov = Rho_Z_InvCov,
                 Rho_Z_CovMu = Rho_Z_CovMu))
}



