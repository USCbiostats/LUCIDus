####### families: record parameters, likelihood #######
normal <- function(K, ...){
  n.par <- 2 * K # number of parameters
  # initialize
  initial.gamma <- function(K, dimCoY){
    structure(list(beta = runif(K + dimCoY),
                   sigma = runif(K)))
  }
  # likelihood function pYgX
  f.pYgX <- function(Y, gamma, K, N, CoY, dimCoY, itr){
    pYgX <- mat.or.vec(N, K)
    mu <- gamma$beta[1:K]
    if(dimCoY != 0 && itr > 1){
      eCoY <- CoY %*% gamma$beta[(K + 1):(K + dimCoY)]
      mu <- sapply(mu, function(x) return(x + eCoY))
    }
    for(i in 1:K){
      if(dimCoY == 0){
        pYgX[, i] <- dnorm(Y, mean = mu[i], sd = gamma$sigma[i])
      } else{
        pYgX[, i] <- sapply(1:N, function(x) return(dnorm(Y[x], mean = mu[x, i], sd = gamma$sigma[i])))
      }
    }
    return(pYgX)
  }
  # update parameters for M step
  f.maxY <- function(Y, r, CoY, K, CoYnames){
    Set0 <- as.data.frame(cbind(Y, r, CoY))
    colnames(Set0) <- c("Y", paste0("LC", 1:K), CoYnames)
    Yfit <- glm(as.formula(paste("Y ~ -1 +", paste(colnames(Set0)[-1], collapse = " + "))), data = Set0, family = gaussian)
    mu <- summary(Yfit)$coefficients[1:K, 1]
    sigma <- sqrt(colSums(r * apply(matrix(mu), 1, function(x){(x - Y)^2})) / colSums(r))
    return(structure(list(beta = mu,
                          sigma = sigma)))
  }
  # switch function, rearrange the parameters
  f.switch <- function(beta, mu, sigma, gamma, K){
    var <- vector(mode = "list", length = K)
    index <- order(gamma$beta[1:K])
    beta <- t(t(beta) - beta[index == 1, ])[index, ]
    mu <- mu[index, ]
    for (i in 1:K) {
      var[[i]] <- sigma[, , i]
    }
    gamma$beta[1:K] <- gamma$beta[index]
    gamma$sigma <- gamma$sigma[index]
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


binary <- function(K, ...){
  n.par <- K # number of parameters
  # initialization EM
  initial.gamma <- function(K){
    matrix(runif(K * 2), nrow = 1)
  }
  # likelihood function pYgX
  f.pYgX <- function(Y, gamma, K, N){
    pYgX <- mat.or.vec(N, K)
    for(i in 1:K){
      p <- exp(gamma[i]) / (1 + exp(gamma[i]))
      pYgX[, i] <- p^Y * (1 - p)^(1 - Y)
    }
    return(pYgX)
  }
  # update parameters for M step
  f.maxY <- function(r, Y, CoY = NULL, K){
    Set0 <- as.data.frame(cbind(Y, r[,-1], CoY))
    colnames(Set0)[2:K] <- paste0("LC", 2:K)
    Yfit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1],collapse="+"))),data=Set0,family="binomial")
    new_gamma[2:K] <- exp(coef(Yfit)[2:K])
  }
  return(structure(list(n.par = n.par,
                        initial.gamma = initial.gamma,
                        f.pYgX = f.pYgX,
                        f.maxY = f.maxY)))
}


poisson <- function(K, ...){
  
}