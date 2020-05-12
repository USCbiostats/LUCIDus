#' @title  Estimate latent unknown clusters with multi-omics data
#' 
#' @description This function estimates the latent clusters by integrating genetic features/environmental exposures, biomarkers with/without the outcome of interest. Variable selection is available for analyzing the high-dimensional data.
#'
#' @param G Genetic features/environmental exposures, a \code{\link{matrix}}.
#' @param Z Biomarkers/other omics data, a \code{\link{matrix}}.
#' @param Y Disease outcome, it is suggested to transform it into a n by 1 \code{\link{matrix}}.
#' @param CoG Optional, matrix. Covariates to be adjusted for estimating the latent cluster.
#' @param CoY Optional, matrix. Covariates to be adjusted for estimating the outcome.
#' @param K Number of latent clusters.
#' @param family Type of outcome Y. It should be choose from "normal", "binary".
#' @param useY Whether or not to include the information of Y to estimate the latent clusters. Default is TRUE.
#' @param control A list of tolerance parameters used by EM algorithm. See \code{\link{def.control}}.
#' @param tune A list of tuning parameters used by variable selection procedure. See \code{\link{def.tune}}
#' @param Z.var.str The variance-covariance structure for the biomarkers. See \code{\link{mclustModelNames}} for details.
#'
#' @return A list which contains the several features of LUCID, including:
#' \item{pars}{Estimates of parameters of LUCID, including beta (estimates of genetic feature/environmental exposure), mu (estimates of cluster-specific biomarker means), sigma (estimates of the cluster-specific biomarker variance-covariance matrix) and gamma(estimates of cluster-specific effect and covariates effect related to the outcome)}
#' \item{K}{Number of latent cluster}
#' \item{Z.var.str}{The model used to estimate the cluster-specific variance-covariance matrix, for further details, see \code{\link{mclust}}}
#' \item{likelihood}{The log likelihood of the LUCID model}
#' \item{post.p}{Predicted probability of belonging to each latent cluster}
#' @importFrom nnet multinom
#' @import mclust
#' @importFrom glmnet glmnet
#' @importFrom glasso glasso
#' @importFrom lbfgs lbfgs
#' @import stats
#' @import utils
#' @export
#' @author Yinqi Zhao, Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Cheng Peng, Jun Wang, Isaac Asante, Stan Louie, Ran Jin, Lida Chatzi, Graham Casey, Duncan C Thomas, David V Conti, A Latent Unknown Clustering Integrating Multi-Omics Data (LUCID) with Phenotypic Traits, Bioinformatics, , btz667, https://doi.org/10.1093/bioinformatics/btz667.
#' @examples
#' \dontrun{
#' set.seed(10)
#' fit1 <- est.lucid(G = G1,Z = Z1,Y = Y1, CoY = CovY, K = 2, family = "binary")
#' fit2 <- est.lucid(G = G1,Z = Z1,Y = Y1, CoY = CovY, K = 2, family = "binary", 
#'   tune = def.tune(Select_Z = TRUE, Rho_Z_InvCov = 0.1, Rho_Z_CovMu = 90, 
#'   Select_G = TRUE, Rho_G = 0.02))
#' }
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
                              model.name = model.best, K = K, ind.na = ind.NA, mu = res.mu)
      if(is.null(new.mu.sigma$mu)){
        print("variable selection failed, restart lucid")
        break
      }
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
        res.mu <- new.mu.sigma$mu
        res.sigma <- new.mu.sigma$sigma
        if(useY){
          res.gamma <- new.gamma
        }
        new.loglik <- sum(log(rowSums(new.likelihood)))
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
  
  res.loglik <- sum(log(rowSums(new.likelihood)))
  pars <- switch_Y(beta = res.beta, mu = res.mu, sigma = res.sigma, gamma = res.gamma, K = K)
  res.r <- res.r[, pars$index]
  colnames(pars$beta) <- c("intercept", Gnames)
  colnames(pars$mu) <- Znames
  if(tune$Select_G == TRUE){
    tt1 <- apply(pars$beta[, -1], 2, range)
    selectG <- abs(tt1[2, ] - tt1[1, ]) > 0.001
  } else{
    selectG <- rep(TRUE, dimG)
  }
  if(tune$Select_Z == TRUE){
    tt2 <- apply(pars$mu, 2, range)
    selectZ <- abs(tt2[2, ] - tt2[1, ]) > 0.001
  } else{
    selectZ <- rep(TRUE, dimZ)
  }
  results <- list(pars = list(beta = pars$beta, mu = pars$mu, sigma = pars$sigma, gamma = pars$gamma),
                  K = K, var.names =list(Gnames = Gnames, Znames = Znames, Ynames = Ynames), Z.var.str = model.best, likelihood = res.loglik, post.p = res.r, family = family,
                  par.control = control, par.tune = tune, select = list(selectG = selectG, selectZ = selectZ), useY = useY)
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
                  G, Z, Y = NULL, family.list, K, N, useY, ind.na, ...){
  pXgG <- pZgX <- pYgX <- matrix(rep(1, N * K), nrow = N)
  if(!is.null(beta)){
    xb <- exp(cbind(rep(1, N), G) %*% t(beta))
    pXgG <- xb/rowSums(xb)
  }
  if(!is.null(mu)){
    for (i in 1:K) {
      pZgX[ind.na != 3, i] <- dmvnorm(Z[ind.na != 3, ], mu[i,], round(sigma[, , i], 9))
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


Mstep_Z <- function(Z, r, selectZ, penalty.mu, penalty.cov,
                    model.name, K, ind.na, mu){
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
        print(paste("glasso failed, restart lucid"))
        break
      }
      else{
        new_sigma[, , k] <- l_cov$w
        new_sigma_inv <- l_cov$wi
        new_sigma_est <- l_cov$w
        try_optim_mu <- try(lbfgs(call_eval = fn, call_grad = gr,
                                  mat = dz, mat2 = dr, k = k,  cov_inv = new_sigma_inv, cov = new_sigma_est,
                                  vars = rep(0, Q), invisible = 1, orthantwise_c = penalty.mu))
        if("try-error" %in% class(try_optim_mu)){
          break
        }
        else{
          new_mu[k, ] <- new_sigma[, , k] %*% (try_optim_mu$par)
        }
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
# use lbfgs to estimate mu with L1 penalty
fn <- function(a, mat, mat2, cov_inv, cov, k){
  Mu <- cov %*% a
  tar <- sum(mat2[, k] * apply(mat, 1, function(v) return(t(v - Mu) %*% cov_inv %*% (v - Mu))))
  return(tar)
}

gr <- function(a, mat, mat2, cov_inv, cov, k){
  Mu <- cov %*% a
  return(2 * apply(mat2[, k] * t(apply(mat, 1, function(v) return(Mu - v))), 2, sum))
}




#' Print the output of \code{est.lucid}
#'
#' @param x An object of LUCID model, returned by \code{est.lucid}
#' @param ... Other arguments to be passed to \code{print}
#' @export
#'
print.lucid <- function(x, ...){
  cat("An object estimated by LUCID model", "\n")
  cat("Outcome type:", x$family, "\n")
  cat("Number of clusters:", "K =", x$K, "\n")
  cat("Variance-Covariance structure for biomarkers:", x$Z.var.str, "model")
}