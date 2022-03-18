binary <- function(K, ...){
  n.par <- K # number of parameters
  # initialize EM
  initial.gamma <- function(K, dimCoY){
    structure(list(beta = runif(K + dimCoY),
                   sigma = NULL))
  }
  # likelihood function pYgX
  f.pYgX <- function(Y, gamma, K, N, CoY, dimCoY, itr){
    pYgX <- mat.or.vec(N, K)
    xb <- gamma$beta[1:K]
    if(dimCoY != 0 && itr > 1){
      eCoY <- CoY %*% gamma$beta[(K + 1):(K + dimCoY)]
      xb <- sapply(xb, function(x) return(x + eCoY))
    }
    p <- exp(xb) / (1 + exp(xb))
    for(i in 1:K){
      if(dimCoY == 0 || itr == 1){
        pYgX[, i] <- dbinom(Y, 1, prob = p[i], log = TRUE)
      } else{
        pYgX[, i] <- dbinom(Y, 1, prob = p[, i], log = TRUE)
      }
    }
    return(pYgX)
  }
  # update parameters for M step
  f.maxY <- function(Y, r, CoY, K, CoYnames){
    if(is.null(CoY)) {
      beta <- apply(r, 2, function(x) return(log(sum(x * Y) / (sum(x) - sum(x * Y)))))
    } else {
      Set0 <- as.data.frame(cbind(Y, r[, -1], CoY))
      colnames(Set0) <- c("Y", paste0("LC", 2:K), CoYnames)
      Yfit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family ="binomial")
      beta <- coef(Yfit) 
      beta[2:K] <- beta[2:K] + beta[1] # log odds for all latent cluster
    }
    return(structure(list(beta = beta,
                          sigma = NULL)))
  }
  # switch function, rearrange the parameters
  f.switch <- function(beta, mu, sigma, gamma, K){
    var <- vector(mode = "list", length = K)
    index <- order(gamma$beta[1:K])
    beta <- t(t(beta) - beta[index[1], ])[index, ]
    mu <- mu[index, ]
    for (i in 1:K) {
      var[[i]] <- sigma[, , i]
    }
    gamma$beta[1:K] <- gamma$beta[1:K][index]
    # transfer log odds back to log OR for LC2 to LCK
    gamma$beta[2:K] <- gamma$beta[2:K] - gamma$beta[1]
    names(gamma$beta)[1:K] <- c("LC1(reference)", paste0("LC", 2:K))
    return(structure(list(beta = beta, 
                          mu = mu,
                          sigma = var,
                          gamma = gamma,
                          index = index)))
  }
  return(structure(list(n.par = n.par,
                        initial.gamma = initial.gamma,
                        f.pYgX = f.pYgX,
                        f.maxY = f.maxY,
                        f.switch = f.switch)))
}


