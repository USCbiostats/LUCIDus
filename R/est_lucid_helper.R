#' Calculate the log-likelihood of cluster assignment for each observation
#'
#' @param beta 
#' @param mu 
#' @param sigma 
#' @param gamma 
#' @param G 
#' @param Z 
#' @param Y 
#' @param family.list 
#' @param K 
#' @param N 
#' @param useY 
#' @param ind.na 
#'
#' @return
#'
Estep <- function(beta, 
                  mu, 
                  sigma, 
                  gamma = NULL,
                  G, 
                  Z, 
                  Y = NULL, 
                  family.list, 
                  K, 
                  N, 
                  useY, 
                  ind.na, ...) {
  # initialize vectors for storing likelihood
  pXgG <- pZgX <- pYgX <- matrix(rep(0, N * K), nrow = N)
  
  # log-likelihood for G -> X
  xb <- cbind(rep(1, N), G) %*% t(beta)
  xb_lse <- apply(xb, 1, lse)
  pXgG <- xb - xb_lse
  
  # log-likelihood for X -> Z
  for (i in 1:K) {
    pZgX[ind.na != 3, i] <- mclust::dmvnorm(data = Z[ind.na != 3, ], 
                                            mean = mu[i,], 
                                            sigma = round(sigma[, , i], 9), 
                                            log = TRUE)
  }
  
  # log-likelihood for X->Y
  if(useY){
    pYgX <- family.list$f.pYgX(Y, gamma, K = K, N = N, ...)
  }
  
  loglik <- pXgG + pZgX + pYgX
  return (loglik)
}





#' Estiamte the association between exposure and latent cluster
#'
#' @param G 
#' @param r 
#' @param selectG 
#' @param penalty 
#' @param dimG 
#' @param K 
#'
Mstep_G <- function(G, r, selectG, penalty, dimG, dimCoG, K) {
  new.beta <- matrix(rep(0, K * (dimG + dimCoG + 1)), nrow = K)
  if(selectG){
    tryLasso <- try(glmnet(as.matrix(G), 
                           as.matrix(r), 
                           family = "multinomial", 
                           lambda = penalty))
    if("try-error" %in% class(tryLasso)){
      breakdown <- TRUE
      print(paste("lasso failed"))
    }
    else{
      new.beta[, 1] <- tryLasso$a0
      new.beta[, -1] <- t(matrix(unlist(lapply(tryLasso$beta, function(x) return(x[,1]))), ncol = K))
    }
  }
  else{
    beta.multilogit <- multinom(as.matrix(r) ~ as.matrix(G))
    new.beta[-1, ] <- coef(beta.multilogit)
  }
  return(new.beta)
}


#' Estimate the mean and variance-covariance structure of GMM part of LUCID
#'
#' @param Z 
#' @param r 
#' @param selectZ 
#' @param penalty.mu 
#' @param penalty.cov 
#' @param model.name 
#' @param K 
#' @param ind.na 
#' @param mu 
#'
Mstep_Z <- function(Z, r, selectZ, penalty.mu, penalty.cov,
                    model.name, K, ind.na, mu) {
  dz <- Z[ind.na != 3, ]
  dr <- r[ind.na != 3, ]
  Q <- ncol(Z)
  new_sigma <- array(rep(0, Q^2 * K), dim = c(Q, Q, K))
  new_mu <- matrix(rep(0, Q * K), nrow = K)
  if(selectZ) {
    k <- 1
    while(k <= K){
      #estimate E(S_k) to be used by glasso
      Z_mu <- t(t(dz) - mu[k, ])
      E_S <- (matrix(colSums(dr[, k] * t(apply(Z_mu, 1, function(x) return(x %*% t(x))))), Q, Q)) / sum(dr[, k])
      #use glasso and E(S_k) to estimate new_sigma and new_sigma_inv
      l_cov <- try(glasso(E_S, penalty.cov))
      if("try-error" %in% class(l_cov)){
        print(paste("glasso failed, restart lucid \n"))
        break
      }
      else{
        new_sigma[, , k] <- l_cov$w
        # function to calculate mean
        new_mu[k, ] <- est.mu(j = k, rho = penalty.mu, z = dz, r = dr, mu = mu[k, ], wi = l_cov$wi)
      }
      k <- k + 1
    }
    if("try-error" %in% class(l_cov)){
      return(structure(list(mu = NULL,
                            sigma = NULL)))
    } else{
      return(structure(list(mu = new_mu,
                            sigma = new_sigma)))
    }
  }
  else{
    z.fit <- mstep(modelName = model.name, data = dz, z = dr)
    return(structure(list(mu = t(z.fit$parameters$mean),
                          sigma = z.fit$parameters$variance$sigma)))
  }
}


#' calculate the log-sum-exp
#' @param vec  a vector of length K
lse <- function(vec) {
  c <- max(vec)
  return(c + log(sum(exp(vec - c))))
}

#' use the log-sum-exp trick to normalize a vector
#' @param vec  a vector of length K
lse_vec <- function(vec) {
  norm_vec <- exp(vec - lse(vec))
  return(norm_vec)
}


#' estimate the penalized mean in M-step Z
#'
#' @param j 
#' @param rho 
#' @param z 
#' @param r 
#' @param mu 
#' @param wi 
#'
est.mu <- function(j, rho, z, r, mu, wi){
  p <- ncol(z)
  res.mu <- rep(0, p)
  mu1 <- sapply(1:p, function(x){
    q1 <- t(t(z) - mu) %*% wi[x, ]
    q2 <- q1 + wi[x, x] * z[, x] - wi[x, x] * (z[, x] - mu[x])
    return(abs(sum(q2 * r[, j])) <= rho)
  })
  mu2 <- sapply(1:p, function(x){
    a <- sum(r[, j] * rowSums(t(wi[x, ] * t(z))))
    b <- sum(r[, j]) * (sum(wi[x, ] * mu) - wi[x, x] * mu[x])
    t1 <- (a - b + rho) / (sum(r[, j]) * wi[x, x]) # mu < 0
    t2 <- (a - b - rho) / (sum(r[, j]) * wi[x, x]) # mu > 0
    if(t1 < 0){
      res <- t1
    } else{
      res <- t2
    }
    return(res)
  })
  res.mu[!mu1] <- mu2[!mu1]
  return(res.mu)
}


#' Print the output of \code{est.lucid}
#'
#' @param x An object of LUCID model, returned by \code{est.lucid}
#' @param ... Other arguments to be passed to \code{print}
#' 
print.lucid <- function(x, ...){
  cat("An object estimated by LUCID model", "\n")
  cat("Outcome type:", x$family, "\n")
  cat("Number of clusters:", "K =", x$K, "\n")
  cat("Variance-Covariance structure for biomarkers:", x$modelName, "model")
}