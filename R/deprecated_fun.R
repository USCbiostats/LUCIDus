# This script documents deprecated functions. All the functions in this ducument
# should be deleted after release of next major update

#' @title Deprecated function est.lucid
#' @description This function deprecates. Please use est_lucid instead.
#' @param G Exposures, a numeric vector, matrix, or data frame. Categorical variable 
#' should be transformed into dummy variables. If a matrix or data frame, rows 
#' represent observations and columns correspond to variables.
#' @param Z Omics data, a numeric matrix or data frame. Rows correspond to observations
#' and columns correspond to variables.
#' @param Y Outcome, a numeric vector. Categorical variable is not allowed. Binary 
#' outcome should be coded as 0 and 1.
#' @param CoG Optional, covariates to be adjusted for estimating the latent cluster.
#' A numeric vector, matrix or data frame. Categorical variable should be transformed 
#' into dummy variables. 
#' @param CoY Optional, covariates to be adjusted for estimating the association 
#' between latent cluster and the outcome. A numeric vector, matrix or data frame. 
#' Categorical variable should be transformed into dummy variables.
#' @param K Number of latent clusters. An integer greater or equal to 2. User 
#' can use \code{\link{lucid}} to determine the optimal number of latent clusters.
#' @param family Distribution of outcome. For continuous outcome, use "normal"; 
#' for binary outcome, use "binary". Default is "normal".
#' @param useY Flag to include information of outcome when estimating the latent 
#' cluster. Default is TRUE.
#' @param tol Tolerance for convergence of EM algorithm. Default is 1e-3.
#' @param max_itr Max number of iterations for EM algorithm.
#' @param max_tot.itr Max number of total iterations for \code{est_lucid} function.
#' \code{est_lucid} may conduct EM algorithm for multiple times if the algorithm 
#' fails to converge.
#' @param Rho_G A scalar. This parameter is the LASSO penalty to regularize
#' exposures. If user wants to tune the penalty, use the wrapper 
#' function \code{lucid}
#' @param Rho_Z_Mu A scalar. This parameter is the LASSO penalty to 
#' regularize cluster-specific means for omics data (Z). If user wants to tune the 
#' penalty, use the wrapper function \code{lucid}
#' @param Rho_Z_Cov A scalar. This parameter is the graphical LASSO
#' penalty to estimate sparse cluster-specific variance-covariance matrices for omics 
#' data (Z). If user wants to tune the penalty, use the wrapper function \code{lucid}
#' @param modelName The variance-covariance structure for omics data. 
#' See \code{mclust::mclustModelNames} for details.
#' @param seed An integer to initialize the EM algorithm or imputing missing values.
#'Default is 123.
#' @param init_impute Method to initialize the imputation of missing values in 
#' LUCID. "mclust" will use \code{mclust:imputeData} to implement EM Algorithm 
#' for Unrestricted General Location Model to impute the missing values in omics 
#' data; \code{lod} will initialize the imputation via relacing missing values by
#' LOD / sqrt(2). LOD is determined by the minimum of each variable in omics data.
#' @param init_par Method to initialize the EM algorithm. "mclust" will use mclust
#' model to initialize parameters; "random" initialize parameters from uniform 
#' distribution.
#' @param verbose A flag indicates whether detailed information for each iteration
#' of EM algorithm is printed in console. Default is FALSE.
#'
#' @export
est.lucid <- function(G, 
                      Z, 
                      Y, 
                      CoG = NULL, 
                      CoY = NULL, 
                      K = 2, 
                      family = c("normal", "binary"), 
                      useY = TRUE, 
                      tol = 1e-3,
                      max_itr = 1e3,
                      max_tot.itr = 1e4,
                      Rho_G = 0,
                      Rho_Z_Mu = 0,
                      Rho_Z_Cov = 0, 
                      modelName = "VVV",
                      seed = 123,
                      init_impute = c("mclust", "lod"),
                      init_par = c("mclust", "random"),
                      verbose = FALSE) {
 stop("function `est.lucid` deprecated. Please use `est_lucid`") 
}


#' @title Deprecated function boot.lucid
#' @description This function deprecates. Please use boot_lucid instead.
#' 
#' @param G Exposures, a numeric vector, matrix, or data frame. Categorical variable 
#' should be transformed into dummy variables. If a matrix or data frame, rows 
#' represent observations and columns correspond to variables.
#' @param Z Omics data, a numeric matrix or data frame. Rows correspond to observations
#' and columns correspond to variables.
#' @param Y Outcome, a numeric vector. Categorical variable is not allowed. Binary 
#' outcome should be coded as 0 and 1.
#' @param CoG Optional, covariates to be adjusted for estimating the latent cluster.
#' A numeric vector, matrix or data frame. Categorical variable should be transformed 
#' into dummy variables. 
#' @param CoY Optional, covariates to be adjusted for estimating the association 
#' between latent cluster and the outcome. A numeric vector, matrix or data frame. 
#' Categorical variable should be transformed into dummy variables.
#' @param model A LUCID model fitted by \code{est.lucid}.
#' @param conf A numeric scalar between 0 and 1 to specify confidence level(s) 
#' of the required interval(s).
#' @param R An integer to specify number of bootstrap replicates for LUCID model.
#' If feasible, it is recommended to set R >= 1000. 
#' @export
boot.lucid <- function(G, 
                       Z, 
                       Y, 
                       CoG = NULL, 
                       CoY = NULL, 
                       model, 
                       conf = 0.95,
                       R = 100) {
  stop("function `boot.lucid` deprecated. Please use `boot_lucid`") 
}


