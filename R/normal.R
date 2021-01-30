# function to estimate normal outcome
normal <- function(K, ...){
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
  f.maxY <- function(Y, r, CoY, K, CoYnames, Tr, beta, sigma){
    Set0 <- as.data.frame(cbind(Y, r[, -1], CoY, Tr))
    Trnames <- colnames(Tr)
    colnames(Set0) <- c("Y", paste0("LC", 2:K), CoYnames, Trnames)
    N <- nrow(r)
    dimCoY <- ncol(CoY)
    if(is.null(Tr)){
      # update gamma
      for(j in 1:K){
        xx <- cbind(matrix(data = rep(0, N * K), nrow = N), CoY)
        xx[, j] <- 1
        sigma[j] <- sqrt(sum((Y - xx %*% beta)^2 * r[, j]) / sum(r[, j]))
      }
      # update beta
      if(is.null(CoY)){
        beta <- sapply(1:K, function(x){sum(r[, x] * Y) / sum(r[, x]) })
      } else{
        for(j in 1:K){
          beta[j] <- sum((Y - CoY %*% beta[(K + 1):(K + dimCoY)]) * r[, j]) / sum(r[, j])
        }
        beta_CoY <- beta[-(1:K)]
        for(l in 1:dimCoY){
          numerator <- sum(sapply(1:K, function(j){
            a <- prod(sigma[-j])^2
            b <- sum(CoY[, l] * (Y - beta[j] - CoY[, -l] %*% beta_CoY[-l]))
            return(a * b)
          }))
          denominator <- sum(sapply(1:K, function(j){
            a <- prod(sigma[-j])^2
            b <- sum(r[, j] * CoY[, l]^2)
            return(a * b)
          }))
          beta[l + K] <- numerator / denominator / 2 ##### why the factor 2?
        }
      }
    } else {
      ################ this part needs to update ##################
      inter <- paste(colnames(Set0)[2:K], ":", Trnames, collapse = " + ")
      cluster <- colnames(Set0)[2:K]
      Yfit <- glm(as.formula(paste("Y ~", paste(c(cluster, Trnames, inter, CoYnames), collapse = " + "))), data = Set0, family = gaussian)
      beta <- summary(Yfit)$coefficients[, 1]
      beta[2:K] <- beta[1] + beta[2:K]
      sigma <- rep(sd(residuals(Yfit)), K)
      ##############################################################
     }
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
  return(structure(list(initial.gamma = initial.gamma,
                        f.pYgX = f.pYgX,
                        f.maxY = f.maxY,
                        f.switch = f.switch)))
}

