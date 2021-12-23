#' Simulate latent cluster X based on exposure
#'
#' @param G an N x M1 matrix represents exposures
#' @param beta an M1 x K matrix represents the effect size between exposures 
#' and latent clusters.
#'
sim_X <-  function(G, beta){
  G <- G - colMeans(G)
  xb <-  G %*% beta
  p <- exp(xb) / rowSums(exp(xb))
  X <-  sapply(1:nrow(xb), function(i){
    rmultinom(1, 1, p[i, ])
  })
  ref <- 1:ncol(beta)
  cl <- sapply(1:ncol(X), function(i){
    return(ref[as.logical(X[, i])])
  })
  return(cl)
}


#' Simulated omics data Z for LUCID model, based on GMM
#'
#' @param X a vector of length N representing the latent cluster
#' @param mu an M2 x K matrix representing mean of GMM
#' @param sigma a list of length K represents the variance covariance structure
#' of GMM, each matrix's size is M2 x M2
#'
sim_Z <- function(X, mu, sigma) {
  Z <- sapply(X, function(cl){
    return(rmvnorm(n = 1, mean = mu[, cl], sigma = sigma[[cl]]))
  })
  return(t(Z))
}



#' Simulate normally distributed outcome Y based on cluster assignment X
#'
#' @param X a vector of length N representing cluster assignment, value ranges 
#' from 1 to K
#' @param CovX a matrix of size N x L representing covariates for X->Y
#' @param beta a vector of length (K + L) for coef X->Y
#' @param sigma a vector of length K, representing the variance for GMM
#' 
sim_Y_normal <- function(X, CovX = NULL, beta, sigma) {
  X <- model.matrix(~as.factor(X) - 1) # convert X into indicator matrix
  X2 <- cbind(X, CovX)
  mu <- X2 %*% beta
  Y <- sapply(1:nrow(mu), function(x){
    rnorm(1, mean = mu[x, ], sd = sigma[as.logical(X[x, ])])
  })
  return(Y)
}



#' Simulate binary outcome Y based on cluster assignment X
#'
#' @param X a vector of length N representing cluster assignment, value ranges 
#' from 1 to K
#' @param CovX a matrix of size N x L representing covariates for X->Y
#' @param beta a vector of length (K - 1 + L) for coef X->Y
#' 
sim_Y_binary <- function(X, CovX = NULL, beta) {
  X <- model.matrix(~as.factor(X)) # convert X into indicator matrix
  XX <- cbind(X, CovX)
  xb <- XX %*% beta
  prob <- exp(xb) / (1 + exp(xb))
  Y <- sapply(prob, function(p) {
    rbinom(1, 1, p)
  })
  return(Y)
}
