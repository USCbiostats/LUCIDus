####### main function: conducting EM algorithm #######
est.lucid <- function(G, Z, Y, CoG = NULL, CoY = NULL, K = 2, family = "normal", 
                      useY = TRUE, control = def.control(), tune = def.tune(), Z.var.str = NULL){
  #### pre-processing ####
  # check data format
  N <- nrow(Y); dimG <- ncol(G); dimZ <- ncol(Z); 
  dimCoG <- ifelse(is.null(CoG), 0, ncol(CoG)); dimCoY <- ifelse(is.null(CoY), 0, ncol(CoY))
  if(is.null(colnames(G))){
    Gnames <- paste0("G", 1:dimG)
  } else {Gnames <- colnames(G)}
  if(is.null(colnames(Z))){
    Znames <- paste0("Z", 1:dimZ)
  } else {Znames <- colnames(Z)}
  if(is.null(colnames(Y))){
    Ynames <- "outcome"
  } else {Ynames <- colnames(Y)}
  CoGnames <- colnames(CoG); CoYnames <- colnames(CoY)
  G <- cbind(G, CoG)
  if(!(is.matrix(G) && is.matrix(Z) && is.matrix(Y))){
    stop("input data should be in the form of matrix")
  }
  if(!is.null(CoY)){
    if(!is.matrix(CoY)){
      stop("input data should be in the form of matrix")
    }
  }
  family.list <- switch(family, normal = normal(K = K, dimCoY), 
                                binary = binary(K = K, dimCoY))
  Mstep_Y <- family.list$f.maxY
  switch_Y <- family.list$f.switch
  # check missing pattern
  ind.NA <- Ind.NA(Z)
  if(sum(ind.NA == 2) != 0){
    # NA.Z <- which(is.na(Z), arr.ind = TRUE)
    Z <- imputeData(Z) # initialize imputation
    Z[ind.NA == 3, ] <- NA
  }

  
  #### conduct the EM algorithm ####
  tot.itr <- 0; convergence <- FALSE
  while(!convergence && tot.itr <= control$max_tot.itr){
    # initialize EM algorithm 
    cat("initialize the LUCID ...", "\n")
    res.beta <- matrix(data = runif(K * (dimG + dimCoG + 1)), nrow = K) 
    res.beta[1, ] <- 0
    invisible(capture.output(mclust.fit <- Mclust(Z[ind.NA != 3, ], G = K)))
    if(is.null(Z.var.str)){
      model.best <- mclust.fit$modelName
    } else{
      model.best <- Z.var.str
    }
    res.mu <- t(mclust.fit$parameters$mean) 
    res.sigma <- mclust.fit$parameters$variance$sigma 
    res.gamma <- family.list$initial.gamma(K, dimCoY)
    res.loglik <- -Inf
    itr <- 0

    while(!convergence && itr <= control$max_itr){
      itr <- itr + 1
      tot.itr <- tot.itr + 1
      check.gamma <-  TRUE
      
      # E step
      new.likelihood <- Estep(beta = res.beta, mu = res.mu, sigma = res.sigma, gamma = res.gamma,
                              G = G, Z = Z, Y = Y, family.list = family.list, itr = itr, CoY = CoY, N = N, K = K, useY = useY, dimCoY = dimCoY, ind.na = ind.NA)
      res.r <- new.likelihood / rowSums(new.likelihood)
      if(!all(is.finite(res.r))){
        cat("iteration", itr,": failed: invalid r, try another seed", "\n")
        break
      } else{
        cat("iteration", itr,": E-step finished.", "\n")
      }
      
      # I step
      # if(sum(ind.NA == 2) != 0 && itr != 1){
      #   Z <- Istep_Z(Z = Z, r = res.r, est.mu = res.mu, ind.na = ind.NA, all.na = NA.Z)
      # }
      
      # M step
      invisible(capture.output(new.beta <- Mstep_G(G = G, r = res.r, selectG = tune$Select_G, penalty = tune$Rho_G, dimG = dimG, K = K)))
      new.mu.sigma <- Mstep_Z(Z = Z, r = res.r, selectZ = tune$Select_Z, penalty.mu = tune$Rho_Z_CovMu, penalty.cov = tune$Rho_Z_InvCov,
                              model.name = model.best, K = K, ind.na = ind.NA)
      if(useY){
        new.gamma <- Mstep_Y(Y = Y, r = res.r, CoY = CoY, K = K, CoYnames)
        check.gamma <- is.finite(unlist(new.gamma))
      }
      
      # control step
      check.value <- all(is.finite(new.beta), is.finite(unlist(new.mu.sigma)), check.gamma)
      singular <- try(sapply(1:K, function(x) return(solve(new.mu.sigma$sigma[, , x]))))
      check.singular <- "try-error" %in% class(singular)
      if(!check.value || check.singular){
        cat("iteration", itr,": Invalid estimates")
        break
      } else{
        res.beta <- new.beta
        res.mu <- t(new.mu.sigma$mu)
        res.sigma <- new.mu.sigma$sigma
        if(useY){
          res.gamma <- new.gamma
        }
        new.loglik <- sum(log(rowSums(sapply(1:N, function(x) return(res.r[x, ] * new.likelihood[x, ])))))
        cat("iteration", itr,": M-step finished, ", "loglike = ", new.loglik, "\n")
        if(abs(res.loglik - new.loglik) < control$tol){
          convergence <- TRUE
          cat("Success: LUCID converges!", "\n")
        }
        res.loglik <- new.loglik
      }
    }
  }
  
  #### summarize the results ####
  if(!useY){
    res.gamma <- Mstep_Y(Y = Y, r = res.r, CoY = CoY, K = K, CoYnames = CoYnames)
  }
  res.likelihood <- Estep(beta = res.beta, mu = res.mu, sigma = res.sigma, gamma = res.gamma,
                          G = G, Z = Z, Y = Y, family.list = family.list, itr = itr, CoY = CoY, N = N, K = K, dimCoY = dimCoY, useY = useY, ind.na = ind.NA)
  res.r <- new.likelihood / rowSums(new.likelihood)
  
  res.loglik <- sum(log(rowSums(sapply(1:N, function(x) return(res.r[x, ] * res.likelihood[x, ])))))
  pars <- switch_Y(beta = res.beta, mu = res.mu, sigma = res.sigma, gamma = res.gamma, K = K)
  res.r <- res.r[, pars$index]
  colnames(pars$beta) <- c("intercept", Gnames)
  results <- list(pars = list(beta = pars$beta, mu = pars$mu, sigma = pars$sigma, gamma = pars$gamma),
                  K = K, var.names =list(Gnames = Gnames, Znames = Znames, Ynames = Ynames), Z.var.str = model.best, likelihood = res.loglik, post.p = res.r, family = family,
                  par.control = control, par.tune = tune)
  class(results) <- c("lucid")
  return(results)
}


####### check the missing patter #######
Ind.NA <- function(Z){
  n <- nrow(Z)
  m <- ncol(Z)
  num.NA <- rowSums(is.na(Z))
  ind.NA <- sapply(1:n, function(x) {return(ifelse(num.NA[x] == 0, 1, 
                                                         ifelse(num.NA[x] == m, 3, 2)))})
  return(ind.NA)
  # 1 = complete, 2 = sporadic, 3 = listwise
}


####### E step: calculate the likelihood #######
Estep <- function(beta = NULL, mu = NULL, sigma = NULL, gamma = NULL,
                  G, Z, Y, family.list, K, N, useY, ind.na, ...){
  pXgG <- pZgX <- pYgX <- matrix(rep(1, N * K), nrow = N)
  if(!is.null(beta)){
    xb <- exp(cbind(rep(1, N), G) %*% t(beta))
    pXgG <- xb/rowSums(xb)
  }
  if(!is.null(mu)){
    for (i in 1:K) {
      pZgX[ind.na != 3, i] <- dmvnorm(Z[ind.na != 3, ], mu[i,], sigma[, , i])
    }
  }
  if(useY){
    pYgX <- family.list$f.pYgX(Y, gamma, K = K, N = N, ...)
  }
  likelihood <- pXgG * pZgX * pYgX
  return (likelihood)
}

####### I step: impute missing values in Z #######
# Istep_Z <- function(Z, r, est.mu, ind.na, all.na){
#   n <- nrow(Z)
#   m <- dim(Z)
#   zr <- colMeans(r[ind.na != 3, ])
#   impute <- as.vector(zr) %*% as.matrix(est.mu)
#   Z[all.na] <- impute[all.na[, 2]]
#   Z[ind.na == 3, ] <- NA
#   return(Z)
# }

####### M step: update the parameters #######
Mstep_G <- function(G, r, selectG, penalty, dimG, K){
  new.beta <- matrix(rep(0, K * (dimG + 1)), nrow = K)
  if(selectG){
    tryLasso <- try(glmnet(as.matrix(G), as.matrix(r), family = "multinomial", lambda = penalty))
    if("try-error" %in% class(tryLasso)){
      breakdown <- TRUE
      print(paste(tot.itr,": lasso failed"))
    }
    else{
      new.beta[,1] <- tryLasso$a0
      new.beta[,-1] <- t(matrix(unlist(lapply(tryLasso$beta, function(x) return(x[,1]))), ncol = K))
    }
  }
  else{
    beta.multilogit <- multinom(as.matrix(r) ~ as.matrix(G))
    new.beta[-1,] <- coef(beta.multilogit)
  }
  return(new.beta)
}


Mstep_Z <- function(Z, r, selectZ, penalty.mu, penalty.cov,
                    model.name, K, ind.na){
  if(selectZ) {
    k <- 1
    while(k <= K && !breakdown){
      #estimate E(S_k) to be used by glasso
      Z_mu <- t(t(Z)-mu[k,])
      E_S <- (matrix(colSums(r[,k]*t(apply(Z_mu,1,function(x) return(x%*%t(x))))),Q,Q))/sum(r[,k])
      #use glasso and E(S_k) to estimate new_sigma and new_sigma_inv
      try_glasso <- try(l_cov <- glasso(E_S,Rho_Z_InvCov))
      if("try-error" %in% class(try_glasso)){
        breakdown <- TRUE
        print(paste(itr,": glasso failed"))
      }
      else{
        new_sigma[[k]] <- l_cov$w
        new_sigma_inv <- l_cov$wi
        new_sigma_est <- l_cov$w
        #use lbfgs to estimate mu with L1 penalty
        fn <- function(a){
          mat <- Z
          cov_inv <- new_sigma_inv
          cov <- new_sigma_est
          Mu <- cov%*%a
          return(sum(r[,k]*apply(mat,1,function(v) return(t(v-Mu)%*%cov_inv%*%(v-Mu)))))
        }
        
        gr <- function(a){
          mat <- Z
          cov <- new_sigma_est
          Mu <- cov%*%a
          return(2.0*apply(r[,k]*t(apply(mat,1,function(v) return(Mu-v))),2,sum))
        }
        
        try_optim_mu <- try(lbfgs(call_eval=fn,call_grad = gr, vars = rep(0,Q), invisible=1, orthantwise_c = Rho_Z_CovMu))
        
        if("try-error" %in% class(try_optim_mu)){
          breakdown <- TRUE
        }
        else{
          new_mu[k,] <- new_sigma[[k]]%*%(try_optim_mu$par)
        }
      }
      k <- k + 1
    }
  }
  else{
    z.fit <- mstep(modelName = model.name, data = Z[ind.na != 3, ], z = r[ind.na != 3, ])
    return(structure(list(mu = z.fit$parameters$mean,
                          sigma = z.fit$parameters$variance$sigma)))
  }
}

####### print, S3 #######
print.lucid <- function(x){
  cat("An object estimated by LUCID model", "\n")
  cat("Outcome type:", x$family, "\n")
  cat("Number of clusters:", "K =", x$K, "\n")
  cat("Variance-Covariance structure for biomarkers:", x$Z.var.str, "model")
}