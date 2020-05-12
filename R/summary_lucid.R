#' Summarize the results of LUCID model
#'
#' @param object A model fitted by \code{\link{est.lucid}}
#' @param boot.se A object returned by \code{\link{boot.lucid}}, which contains the bootstrap standard error
#' @param ... Other paramters to be passed to \code{summary}
#' @return A list with class "sumlucid", which contains the following object
#' \item{Beta}{Estimates of genetic/environmental effects (and effect of covariates if included), matrix}
#' \item{Mu}{Estimates of cluster-specific biomarker means, matrix}
#' \item{Gamma}{Estimates of cluster-specific disease risk (and effect of covariates if included), vector}
#' \item{Family}{Type of Y, binary or normal}
#' \item{K}{Number of latent clusters}
#' \item{loglik}{log likelihood of the model}
#' \item{BIC}{Bayesian Information Criteria of the model}
#' \item{boot.se}{Bootstrap SE for estimates, an object returned by \code{\link{boot.lucid}}}
#' @export
#' @author Yinqi Zhao, Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Cheng Peng, Jun Wang, Isaac Asante, Stan Louie, Ran Jin, Lida Chatzi, Graham Casey, Duncan C Thomas, David V Conti, A Latent Unknown Clustering Integrating Multi-Omics Data (LUCID) with Phenotypic Traits, Bioinformatics, , btz667, https://doi.org/10.1093/bioinformatics/btz667.
#' 
#' @examples
#' \dontrun{
#' fit1 <- est.lucid(G = G1, Z = Z1, Y = Y1, CoY = CovY, K = 2, family = "binary", useY = FALSE)
#' summary(fit1)
#' fit2 <- est.lucid(G = G1, Z = Z1, Y = Y1, CoY = CovY, K = 2, family = "binary", useY = FALSE,
#'  tune = def.tune(Select_Z = TRUE, Rho_Z_InvCov = 0.1, Rho_Z_CovMu = 90, 
#'  Select_G = TRUE, Rho_G = 0.02)) 
#' summary(fit2)
#' }

summary.lucid <- function(object, boot.se = NULL, ...){
  s1 <- object$select$selectG
  s2 <- object$select$selectZ
  nG <- sum(s1)
  nZ <- sum(s2)
  K <- object$K
  gamma <- object$pars$gamma
  if(object$family == "normal"){
    nY <- length(gamma$beta) + length(gamma$sigma)
  }
  if(object$family == "binary"){
    nY <- length(gamma$beta)
  }
  npars <- (nG + 1) * (K - 1) + (nZ * K + nZ * (nZ + 1) / 2 * K) + nY
  BIC <- -2 * object$likelihood + npars * log(nrow(object$post.p))
  results <- list(beta = object$pars$beta[, c(TRUE, s1)],
                    mu = object$pars$mu[, s2],
                    gamma = object$pars$gamma,
                    family = object$family,
                    K = K,
                    BIC = BIC,
                  loglik = object$likelihood,
                  boot.se = boot.se)
  class(results) <- "sumlucid"
  return(results)
}



#' Print the output of LUCID in a nicer table
#'
#' @param x An object returaned by \code{summary.lucid}
#' @param ... Other paramters to be passed to \code{print}
#' @export
#'
print.sumlucid <- function(x, ...){
  K <- x$K
  beta <- x$beta
  dim1 <- ncol(beta) - 1
  z.mean <- as.data.frame(t(x$mu))
  cat("----------Summary of the LUCID model---------- \n \n")
  cat("K = ", K, ", log likelihood =", x$loglik, ", BIC = ", x$BIC, "\n \n")
  y <- switch(x$family, normal = f.normal,
              binary = f.binary)
  y(x$gamma, K, se = x$boot.se$gamma)
  cat("\n")
  cat("(2) Z: estimates of biomarker means for each latent cluster \n")
  if(is.null(x$boot.se)){
    colnames(z.mean) <- paste0("cluster", 1:K)
    print(z.mean)
  } else{
    print(x$boot.se$mu)
  }
  cat("\n")
  cat("(3) E: the odds ratio of being assigned to each latent cluster for each exposure \n")
  dd <- as.matrix(as.data.frame(beta)[2:K, 2:ncol(beta)])
  g.or <- data.frame(original = unlist(split(dd, row(dd))))
  rownames(g.or) <- paste0(colnames(beta)[-1], ".cluster", sapply(2:K, function(x) return(rep(x, dim1))))
  if(is.null(x$boot.se)){
    g.or$OR <- exp(g.or$original)
    print(g.or)
  } else{
    bb <- x$boot.se$beta
    g.or <- cbind(g.or, bb[, 2:4], OR = exp(bb[, 1]), OR.L = exp(bb[, 3]), OR.U = exp(bb[, 4]))
    print(g.or)
  }
}


f.normal <- function(x, K, se){
  gamma <- x$beta
  sigma <- x$sigma
  cat("(1) Y (normal outcome): the mean and the sd of the Gaussian mixture Y for each latent cluster \n")
  y <- matrix(c(gamma, sigma), ncol = 2)
  row.names(y) <- paste0("cluster", 1:K)
  colnames(y) <- c("mu", "sd")
  if(!is.null(se)){
    gamma <- cbind(y, se[, 2:4])
  }
  print(y)
}

f.binary <- function(x, K, se){
  cat("(1) Y (binary outcome): odds ratio of Y for each latent cluster (covariate) \n")
  gamma <- as.data.frame(x$beta)
  colnames(gamma) <- "Original"
  if(is.null(se)){
    gamma$OR <- exp(gamma$Original)
  } else{
    gamma <- cbind(gamma, se[, 2:4], OR = exp(gamma[, 1]), OR.L = exp(se[, 3]), OR.U = exp(se[, 4]))
  }
  print(gamma)
}


