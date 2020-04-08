summary.lucid <- function(x){
  K <- x$K
  beta <- x$pars$beta
  z.mean <- as.data.frame(t(x$pars$mu))
  cat("~Summary of the LUCID model~ \n \n")
  cat("K = ", K, ", log likelihood =", x$likelihood, "\n \n")
  y <- switch(x$family, normal = f.normal,
                        binary = f.binary)
  y(x$pars$gamma, K)
  cat("\n")
  cat("(2) Z: effect of biomarkers for each latent cluster \n")
  colnames(z.mean) <- paste0("cluster", 1:K)
  print(z.mean)
  cat("\n")
  cat("(3) E: the odds ratio of being assigned to each latent cluster for each exposure \n")
  g.or <- t(exp(beta)[, 2:ncol(beta)])
  colnames(g.or) <- paste0("cluster", 1:K)
  print(g.or)
}



f.normal <- function(x, K){
  gamma <- x$beta
  sigma <- x$sigma
  cat("(1) Y: the mean expression and the variance of Y for each latent cluster \n")
  y <- matrix(c(gamma, sigma), ncol = 2)
  row.names(y) <- paste0("cluster", 1:K)
  colnames(y) <- c("mean", "variance")
  print(y)
}

f.binary <- function(x, K){
  gamma <- as.data.frame(x$beta)
  cat("(1) Y (binary outcome): odds ratio of Y for each latent cluster (covariate) \n")
  colnames(gamma) <- "Estimate"
  row.names(gamma)[1:K] <- paste0("cluster", 1:K)
  print(exp(gamma))
}


