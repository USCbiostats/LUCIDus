
normal <- function(K, ...){
  n.par <- 2 * K # number of parameters
  # initialize EM
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
      if(dimCoY == 0 || itr == 1){
        pYgX[, i] <- dnorm(Y, mean = mu[i], sd = gamma$sigma[i])
      } else{
        pYgX[, i] <- sapply(1:N, function(x) return(dnorm(Y[x], mean = mu[x, i], sd = gamma$sigma[i])))
      }
    }
    return(pYgX)
  }
  # update parameters for M step
  f.maxY <- function(Y, r, CoY, K, CoYnames, Tr){

      Set0 <- as.data.frame(cbind(Y, r[, -1], CoY, Tr))
      Trnames <- colnames(Tr)
      colnames(Set0) <- c("Y", paste0("LC", 2:K), CoYnames, Trnames)
      if(is.null(Tr)){
        Yfit <- glm(as.formula(paste("Y ~", paste(colnames(Set0)[-1], collapse = " + "))), data = Set0, family = gaussian)
      } else {
        inter <- paste(colnames(Set0)[2:K], ":", Trnames, collapse = " + ")
        cluster <- colnames(Set0)[2:K]
        Yfit <- glm(as.formula(paste("Y ~", paste(c(cluster, Trnames, inter, CoYnames), collapse = " + "))), data = Set0, family = gaussian)
      }
      beta <- summary(Yfit)$coefficients[, 1]
      beta[2:K] <- beta[1] + beta[2:K]
      sigma <- rep(sd(residuals(Yfit)), K)
    # if(!is.null(CoY)){
    # } else{
    #   beta <- sapply(1:K, function(x){sum(r[, x] * Y) / sum(r[, x]) })
    #   sigma <- sqrt(colSums(r * apply(matrix(beta), 1, function(x){(x - Y)^2})) / colSums(r))
    # }
    return(structure(list(beta = beta,
                          sigma = sigma)))
  }
  # switch function, rearrange the parameters
  f.switch <- function(beta, mu, sigma, gamma, K, Tr){
    var <- vector(mode = "list", length = K)
    index <- order(gamma$beta[1:K])
    beta <- t(t(beta) - beta[index[1], ])[index, ]
    mu <- mu[index, ]
    for (i in 1:K) {
      var[[i]] <- sigma[, , i]
    }
    gamma$beta[1:K] <- gamma$beta[index]
    if(!is.null(Tr)){
      int <- c(0, gamma$beta[(K + 2):(2 * K)])
      gamma$beta[K + 1] <- gamma$beta[K + 1] + int[index == 1]
      int <- (int - int[index == 1])[index]
      gamma$beta[(K + 2):(2 * K)] <- int[-1]
    }
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
        pYgX[, i] <- p[i]^Y * (1 - p[i])^(1 - Y)
      } else{
        pYgX[, i] <- sapply(1:N, function(x) return(p[x, i]^Y[x] * (1 - p[x, i])^(1 - Y[x])))
      }
    }
    return(pYgX)
  }
  # update parameters for M step
  f.maxY <- function(Y, r, CoY, K, CoYnames, Tr){
    Set0 <- as.data.frame(cbind(Y, r[, -1], CoY, Tr))
    Trnames <- colnames(Tr)
    colnames(Set0) <- c("Y", paste0("LC", 2:K), CoYnames, Trnames)
    if(is.null(Tr)){
      Yfit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family = binomial)
    } else {
      inter <- paste(colnames(Set0)[2:K], ":", Trnames, collapse = " + ")
      cluster <- colnames(Set0)[2:K]
      Yfit <- glm(as.formula(paste("Y ~", paste(c(cluster, Trnames, inter, CoYnames), collapse = " + "))), data = Set0, family = binomial)
    }
    beta <- coef(Yfit) # this is the baseline log odds
    beta[2:K] <- beta[1] + beta[2:K] # log odds for each latent cluster
    return(structure(list(beta = beta,
                          sigma = NULL)))
  }
  # switch function, rearrange the parameters
  f.switch <- function(beta, mu, sigma, gamma, K, Tr){
    var <- vector(mode = "list", length = K)
    index <- order(gamma$beta[1:K])
    beta <- t(t(beta) - beta[index[1], ])[index, ]
    mu <- mu[index, ]
    for (i in 1:K) {
      var[[i]] <- sigma[, , i]
    }
    ref <- gamma$beta[1:K][index == 1]
    gamma$beta[1:K] <- (gamma$beta[1:K] - ref)[index]
    gamma$beta[1] <- ref
    if(!is.null(Tr)){
      int <- c(0, gamma$beta[(K + 2):(2 * K)])
      gamma$beta[K + 1] <- gamma$beta[K + 1] + int[index == 1]
      int <- (int - int[index == 1])[index]
      gamma$beta[(K + 2):(2 * K)] <- int[-1]
    }
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



############# need more test: 11/06/20 ####################
poisson <- function(K, ...){
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
    lambda <- exp(xb)
    for(i in 1:K){
      if(dimCoY == 0 || itr == 1){
        pYgX[, i] <- exp(-lambda[i]) * lambda[i]^Y / factorial(Y)
      } else{
        pYgX[, i] <- sapply(1:N, function(x) return(exp(-lambda[x, i]) * lambda[x, i]^Y[x] / factorial(Y[x])))
      }
    }
    return(pYgX)
  }
  # update parameters for M step
  f.maxY <- function(Y, r, CoY, K, CoYnames){
    Set0 <- as.data.frame(cbind(Y, r[, -1], CoY))
    colnames(Set0) <- c("Y", paste0("LC", 2:K), CoYnames)
    Yfit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family = poisson)
    beta <- coef(Yfit) # this is the baseline log odds
    beta[2:K] <- beta[1] + beta[2:K] # log odds for each latent cluster
    # beta <- c(0, coef(Yfit)[-1])
    return(structure(list(beta = beta,
                          sigma = NULL)))
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
    ref <- gamma$beta[1:K][index == 1]
    gamma$beta[1:K] <- (gamma$beta[1:K] - ref)[index]
    gamma$beta[1] <- ref
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