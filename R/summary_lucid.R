summary.lucid <- function(x){
  K <- x$K
  z.mean <- x$pars$mu
  beta <- x$pars$beta
  cat("~Summary of the LUCID model~ \n \n")
  cat("K = ", K, ", log likelihood =", x$likelihood, "\n \n")
  y <- switch(x$family, normal = f.normal,
                        binary = f.binary)
  y(x$pars$gamma)
  cat("\n")
  cat("(2) Z: the mean expression of biomarkers for each latent cluster \n")
  if(is.null(x$var.names$Znames)){
    colnames(z.mean) <- paste0("Z", 1:ncol(z.mean))
  } else{
    colnames(z.mean) <- x$var.names$Znames
  }
  row.names(z.mean) <- paste0("cluster", 1:K)
  print(z.mean)
  cat("\n")
  cat("(3) E: the odds ratio of being assigned to each latent cluster for each exposure \n")
  g.or <- exp(beta)[, 2:ncol(beta)]
  if(is.null(x$var.names$Gnames)){
    colnames(g.or) <- c(paste0("G", 1:ncol(g.or)))
  } else{
    colnames(g.or) <- x$var.names$Gnames
  }
  row.names(g.or) <- c("cluster1(reference)", paste0("cluster", 2:K))
  print(g.or)
}



f.normal <- function(x){
  gamma <- x$beta
  sigma <- x$sigma
  n <- length(gamma)
  cat("(1) Y: the mean expression and the variance of Y for each latent cluster \n")
  y <- matrix(c(gamma, sigma), ncol = 2)
  row.names(y) <- paste0("Cluster", 1:n)
  colnames(y) <- c("mean", "variance")
  print(y)
}

f.binary <- function(x){
  
}


