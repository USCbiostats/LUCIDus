#' @title Summarize results of LUCID model
#'
#' @param object A LUCID model fitted by \code{\link{est.lucid}}
#' @param boot.se An object returned by \code{\link{boot.lucid}}, 
#' which contains the bootstrap confidence intervals
#' 
#' @export
#' @examples
#' \dontrun{
#' # use simulated data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#' 
#' # fit lucid model
#' fit1 <- est.lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, 
#' seed = 1008)
#' 
#' # conduct bootstrap resampling
#' boot1 <- boot.lucid(G = G, Z = Z, Y = Y_normal, model = fit1, R = 100)
#' 
#' # summarize lucid model
#' summary_lucid(fit1)
#' 
#' # summarize lucid model with bootstrap CIs
#' summary_lucid(fit1, boot.se = boot1)
#' }

summary_lucid <- function(object, boot.se = NULL){
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
  npars <- (nG + 1) * (K - 1) + (nZ * K + nZ^2 * K) + nY
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
#' @param x An object returned by \code{summary.lucid}
#' @param ... Other parameters to be passed to \code{print}
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
  cat("(2) Z: mean of omics data for each latent cluster \n")
  if(is.null(x$boot.se)){
    colnames(z.mean) <- paste0("mu_cluster", 1:K)
    print(z.mean)
  } else{
    print(x$boot.se$mu)
  }
  cat("\n")
  cat("(3) E: odds ratio of being assigned to each latent cluster for each exposure \n")
  dd <- as.matrix(as.data.frame(beta)[2:K, 2:ncol(beta)])
  g.or <- data.frame(beta = unlist(split(dd, row(dd))))
  rownames(g.or) <- paste0(colnames(beta)[-1], ".cluster", sapply(2:K, function(x) return(rep(x, dim1))))
  if(is.null(x$boot.se)){
    g.or$OR <- exp(g.or$beta)
    print(g.or)
  } else{
    # bb <- x$boot.se$beta
    # g.or <- cbind(g.or, bb[, 2:4], OR = exp(bb[, 1]), OR.L = exp(bb[, 3]), OR.U = exp(bb[, 4]))
    print(x$boot.se$beta)
  }
}


# summarize output of normal outcome
f.normal <- function(x, K, se){
  
  cat("(1) Y (normal outcome): the mean of Y for each latent cluster (and effect of covariates) \n")
  
  if(!is.null(se)){
    y <- se
  } else {
    gamma <- x$beta
    y <- as.data.frame(gamma)
    row.names(y)[1:K] <- paste0("cluster", 1:K)
    colnames(y) <- "beta"
  }
  print(y)
}


# summarize output of binary outcome
f.binary <- function(x, K, se){
  cat("(1) Y (binary outcome): log odds of Y for each latent cluster (and log OR of covariate)\n")
  gamma <- as.data.frame(x$beta)
  colnames(gamma) <- "gamma"
  if(is.null(se)){
    gamma$`exp(gamma)` <- exp(gamma$gamma)
  } else{
    gamma <- cbind(gamma, se[, -1])
  }
  print(gamma)
}

