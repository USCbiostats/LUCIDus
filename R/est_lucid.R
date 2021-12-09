#' @title  Estimate latent unknown clusters with multi-omics data
#' 
#' @description This function fits a LUCID model, which estimates the latent clusters by integrating genetic features/environmental exposures, biomarkers with/without the outcome of interest. Variable selection is available for analyzing the high-dimensional data.
#'
#' @param G Exposures, a \code{matrix} or a \code{dataframe}.
#' @param Z Omics data, a \code{matrix} or a \code{dataframe}.
#' @param Y Outcome, an n x 1 \code{matrix}.
#' @param CoG Optional, a \code{matrix}. Covariates to be adjusted for estimating the latent cluster.
#' @param CoY Optional, a \code{matrix}. Covariates to be adjusted for estimating the association between latent cluster and outcome.
#' @param K The number of latent clusters. An integer greater or equal to 2. You can use \code{\link{tune.lucid}} to determine the optimal number of latent clusters.
#' @param family Distribution of outcome Y. Currently accepted distributions are "normal" and "binary".
#' @param useY Whether or not to include the information of Y to estimate the latent clusters. Default is TRUE.
#' @param control A list of tolerance parameters used by EM algorithm. See \code{\link{def.control}}.
#' @param tune A list of tuning parameters used by variable selection procedure. See \code{\link{def.tune}}
#' @param modelName The variance-covariance structure for omics data. See \code{\link{mclustModelNames}} for details.
#' @param seed An integer to initialize the EM algorithm
#' 
#' @return A list which contains the several features of LUCID, including:
#' \item{pars}{Estimates of parameters of LUCID, including beta (estimates of genetic feature/environmental exposure), mu (estimates of cluster-specific biomarker means), sigma (estimates of the cluster-specific biomarker variance-covariance matrix) and gamma(estimates of cluster-specific effect and covariates effect related to the outcome)}
#' \item{K}{Number of latent cluster}
#' \item{modelName}{The model used to estimate the cluster-specific variance-covariance matrix, for further details, see \code{\link{mclust}}}
#' \item{likelihood}{The log likelihood of the LUCID model}
#' \item{post.p}{Predicted probability of belonging to each latent cluster}
#' \item{Z}{}
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
#' fit1 <- est.lucid(G = sim1[, 1:10],Z = sim1[, 11:20],Y = as.matrix(sim1[, 21]), 
#' K = 2, family = "binary")
#' fit2 <- est.lucid(G = sim1[, 1:10],Z = sim1[, 11:20],Y = as.matrix(sim1[, 21]), 
#' K = 2, family = "binary", 
#' tune = def.tune(Select_Z = TRUE, Rho_Z_InvCov = 0.1, Rho_Z_CovMu = 90, 
#' Select_G = TRUE, Rho_G = 0.02))
#' }
est.lucid <- function(G, Z, Y, 
                      CoG = NULL, CoY = NULL, 
                      K = 2, 
                      family = c("normal", "binary"), 
                      useY = TRUE, 
                      control = def.control(), 
                      tune = def.tune(), 
                      modelName = NULL,
                      seed = 123) {
  
  #============ 1. basic setup for estimation function =============
  family <- match.arg(family)
  
  # 1.1 check data format 
  if(is.data.frame(G) | is.vector(G)) {
    G <- as.matrix(G)
    if(!is.numeric(G)) {
      stop("'G' must be numeric (for categorical input, use dummy variables)")
    }
  }
  if(is.null(colnames(G))){
    Gnames <- paste0("G", 1:ncol(G))
  } else {Gnames <- colnames(G)}
  
  if(is.data.frame(Z) | is.vector(Z)) {
    Z <- as.matrix(Z)
    if(!is.numeric(Z)) {
      stop("'Z' must be numeric (for categorical input, use dummy variables)")
    }
  }
  if(is.null(colnames(Z))){
    Znames <- paste0("Z", 1:ncol(Z))
  } else {Znames <- colnames(Z)}
  
  if(is.data.frame(Y) | is.vector(Y)) {
    Y <- as.matrix(Y)
    if(!is.numeric(Y)) {
      stop("'Y' must be numeric (for categorical input, use dummy variables)")
    }
  }
  
  if(is.null(colnames(Y))){
    Ynames <- "outcome"
  } else {Ynames <- colnames(Y)}
  
  if(!is.null(CoG)) {
    if(is.data.frame(CoG) | is.vector(CoG)) {
      CoG <- as.matrix(CoG)
      if(!is.numeric(CoG)) {
        stop("'CoG' must be numeric (for categorical input, use dummy variables)")
      }
    }
    CoGnames <- colnames(CoG)
  }
  if(!is.null(CoY)) {
    if(is.data.frame(CoY) | is.vector(CoY)) {
      CoY <- as.matrix(CoY)
      if(!is.numeric(CoG)) {
        stop("'CoY' must be numeric (for categorical input, use dummy variables)")
      }
    }
    CoYnames <- colnames(CoY)
  }
  
  # 1.2 record input dimensions, family function 
  N <- nrow(Y)
  dimG <- ncol(G)
  dimZ <- ncol(Z); 
  dimCoG <- ifelse(is.null(CoG), 0, ncol(CoG))
  dimCoY <- ifelse(is.null(CoY), 0, ncol(CoY))
  G <- cbind(G, CoG)
  family.list <- switch(family, normal = normal(K = K, dimCoY), 
                        binary = binary(K = K, dimCoY))
  Mstep_Y <- family.list$f.maxY
  switch_Y <- family.list$f.switch
  
  
  # 1.3. check missing pattern
  na_pattern <- check_na(Z)
  if(na_pattern$impute_flag) {
    # initialize imputation
    cat("Intializing the missing values in omics data 'Z'\n")
    invisible(capture.output(Z <- mclust::imputeData(Z)))
    Z[na_pattern$indicator_na == 3, ] <- NA
  }
  
  
  #=============== 2. EM algorithm for LUCID ================
  tot.itr <- 0
  convergence <- FALSE
  while(!convergence && tot.itr <= control$max_tot.itr) {
    if(tot.itr > 0) {
      seed <- seed + 10
    }
    set.seed(seed)
    
    # 2.1 initialize model parameters
    cat("initialize LUCID with mclust \n")
    
    # initialize beta 
    res.beta <- matrix(data = runif(K * (dimG + dimCoG + 1)), nrow = K) 
    res.beta[1, ] <- 0
    
    # initialize GMM parameters mu and sigma 
    invisible(capture.output(mclust.fit <- Mclust(Z[na_pattern$indicator_na != 3, ], 
                                                  G = K,
                                                  modelNames = modelName)))
    if(is.null(modelName)){
      model.best <- mclust.fit$modelName
    } else{
      model.best <- modelName
    }
    res.mu <- t(mclust.fit$parameters$mean)
    res.sigma <- mclust.fit$parameters$variance$sigma
    
    # initialize family specific parameters gamma
    res.gamma <- family.list$initial.gamma(K, dimCoY)
    
    
    # start EM algorithm
    res.loglik <- -Inf
    itr <- 0
    while(!convergence && itr <= control$max_itr){
      itr <- itr + 1
      tot.itr <- tot.itr + 1
      check.gamma <-  TRUE
      
      # 2.2 E-step:
      # calculate log-likelihood for observation i being assigned to cluster j
      new.likelihood <- Estep(beta = res.beta, 
                              mu = res.mu, 
                              sigma = res.sigma, 
                              gamma = res.gamma,
                              G = G, Z = Z, Y = Y, CoY = CoY, 
                              N = N, 
                              K = K, 
                              family.list = family.list, 
                              itr = itr, 
                              useY = useY, 
                              dimCoY = dimCoY, 
                              ind.na = na_pattern$indicator_na)
      # normalize the log-likelihood to probability
      res.r <- t(apply(new.likelihood, 1, lse_vec))

      if(!all(is.finite(res.r))){
        cat("iteration", itr,": EM algorithm collapsed: invalid estiamtes due to over/underflow, try another seed \n")
        break
      } else{
        cat("iteration", itr,": E-step finished.\n")
      }
      
      
      # 2.3 M-step - parameters
      # update model parameters to maximize the expected likelihood
      invisible(capture.output(new.beta <- Mstep_G(G = G, 
                                                   r = res.r, 
                                                   selectG = tune$Select_G, 
                                                   penalty = tune$Rho_G, 
                                                   dimG = dimG, 
                                                   K = K)))
      new.mu.sigma <- Mstep_Z(Z = Z, 
                              r = res.r, 
                              selectZ = tune$Select_Z, 
                              penalty.mu = tune$Rho_Z_CovMu, 
                              penalty.cov = tune$Rho_Z_InvCov,
                              model.name = model.best, 
                              K = K, 
                              ind.na = na_pattern$indicator_na, 
                              mu = res.mu)
      if(is.null(new.mu.sigma$mu)){
        print("variable selection failed, restart lucid \n")
        break
      }
      if(useY){
        new.gamma <- Mstep_Y(Y = Y, r = res.r, CoY = CoY, K = K, CoYnames)
        check.gamma <- is.finite(unlist(new.gamma))
      }
      
      
      # 2.4 M step - impute missing values
      if(na_pattern$impute_flag){
        Z <- Istep_Z(Z = Z, 
                     p = res.r, 
                     mu = res.mu, 
                     sigma = res.sigma,
                     index = na_pattern$index)
      }
      
      
      # 2.5 control step
      check.value <- all(is.finite(new.beta), is.finite(unlist(new.mu.sigma)), check.gamma)
      singular <- try(sapply(1:K, function(x) return(solve(new.mu.sigma$sigma[, , x]))))
      check.singular <- "try-error" %in% class(singular)
      if(!check.value || check.singular){
        cat("iteration", itr,": Invalid estimates \n")
        break
      } else{
        res.beta <- new.beta
        res.mu <- new.mu.sigma$mu
        res.sigma <- new.mu.sigma$sigma
        if(useY){
          res.gamma <- new.gamma
        }

        new.loglik <- sum(log(rowSums(exp(new.likelihood))))

        if(tune$Select_G) {
          new.loglik <- new.loglik - tune$Rho_G * sum(abs(res.beta))
        }
        if(tune$Select_Z) {
          new.loglik <- new.loglik - tune$Rho_Z_CovMu * sum(abs(res.mu)) - tune$Rho_Z_InvCov * sum(abs(res.sigma))
        }
        if(tune$Select_G | tune$Select_Z) {
          cat("iteration", itr,": M-step finished, ", "penalized loglike = ", sprintf("%.3f", new.loglik), "\n")
        } else{
          cat("iteration", itr,": M-step finished, ", "loglike = ", sprintf("%.3f", new.loglik), "\n")
        }
        if(abs(res.loglik - new.loglik) < control$tol){
          convergence <- TRUE
          cat("Success: LUCID converges!", "\n")
        }
        res.loglik <- new.loglik
      }
    }
  }
  
  #============ 3. summarize results ===============
  if(!useY){
    res.gamma <- Mstep_Y(Y = Y, r = res.r, CoY = CoY, K = K, CoYnames = CoYnames)
  }
  res.likelihood <- Estep(beta = res.beta, 
                          mu = res.mu, 
                          sigma = res.sigma, 
                          gamma = res.gamma,
                          G = G, Z = Z, Y = Y, 
                          family.list = family.list, 
                          itr = itr, CoY = CoY, N = N, K = K, 
                          dimCoY = dimCoY, useY = useY, 
                          ind.na = na_pattern$indicator_na)
  res.r <- t(apply(res.likelihood, 1, lse_vec))
  

  res.loglik <- sum(log(rowSums(exp(res.likelihood))))

  if(tune$Select_G) {
    res.loglik <- res.loglik - tune$Rho_G * sum(abs(res.beta))
  }
  if(tune$Select_Z) {
    res.loglik <- res.loglik - tune$Rho_Z_CovMu * sum(abs(res.mu)) - tune$Rho_Z_InvCov * sum(abs(res.sigma))
  }
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
  results <- list(pars = list(beta = pars$beta, 
                              mu = pars$mu, 
                              sigma = pars$sigma, 
                              gamma = pars$gamma),
                  K = K, 
                  var.names =list(Gnames = Gnames, 
                                  Znames = Znames, 
                                  Ynames = Ynames), 
                  modelName = model.best, 
                  likelihood = res.loglik,
                  post.p = res.r, 
                  family = family,
                  par.control = control, 
                  par.tune = tune, 
                  select = list(selectG = selectG, selectZ = selectZ), 
                  useY = useY,
                  Z = Z)
  class(results) <- c("lucid")
  return(results)
}



#' Check missing patterns
#'
#' @param Z 
#'
#' @return
#' @examples
check_na <- function(Z){
  N <- nrow(Z)
  M <- ncol(Z)
  index <- !is.na(Z)
  obs_na <- rowSums(!index)
  indicator_na <- sapply(1:N, function(i) {
    return(ifelse(obs_na[i] == 0, 1,
                  ifelse(obs_na[i] == M, 3, 2)))
  })
  impute_flag <- sum(indicator_na == 2) != 0
  # 1 = complete, 2 = sporadic, 3 = listwise

  return(list(index = index,
              indicator_na = indicator_na,
              impute_flag = impute_flag))
}





#' calculate the log-likelihood of cluster assignment for each observation
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
#' @export
#'
#' @examples
Estep <- function(beta, mu, sigma, gamma,
                  G, Z, Y = NULL, family.list, K, N, useY, ind.na, ...) {
  pXgG <- pZgX <- pYgX <- matrix(rep(1, N * K), nrow = N)
  
  # log-likelihood for G->X
  xb <- cbind(rep(1, N), G) %*% t(beta)
  pXgG <- log(t(apply(xb, 1, lse_vec)))
  
  # log-likelihood for X->Z
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


#' use the log-sum-exponential trick to avoid over/underflow
#' @param vec  a vector of length K
lse_vec <- function(vec) {
  c <- max(vec)
  lse <- log(sum(exp(vec - c)))
  return(exp(vec - c - lse))
}




#' Title
#'
#' @param G 
#' @param r 
#' @param selectG 
#' @param penalty 
#' @param dimG 
#' @param K 
#'
#' @return
#' @export
#'
#' @examples
Mstep_G <- function(G, r, selectG, penalty, dimG, K) {
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


#' Title
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
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param j 
#' @param rho 
#' @param z 
#' @param r 
#' @param mu 
#' @param wi 
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param Z 
#' @param r 
#' @param mu
#' @param sigma
#' @param index
#'
#' @return
#' @export
#'
#' @examples
Istep_Z <- function(Z, p, mu, sigma, index){
  N <- nrow(Z)
  Z_fill <- t(sapply(1:N, function(i) {
    fill_data(obs = Z[i, ], mu = mu, sigma = sigma, p = p, index = index[i, ])
  }))
  return(Z_fill)
}


#' Fill in missing data by optimizing the likelihood function
#'
#' @param obs a vector of length M
#' @param mu a matrix of size M x K
#' @param sigma a matrix of size M x M x K
#' @param p a vector of length K
#' @param index a vector of length M, indicating whether a value is missing 
#' or not in the raw data
#'
#' @return an observation with updated imputed value
#'
fill_data <- function(obs, mu, sigma, p, index) {
  mu <- t(mu)
  M <- length(obs)
  K <- ncol(mu)
  # impute missing values
  if(any(!index) && !all(!index)) {
    sigma_inv <- array(rep(0, M * M * K), dim = c(M, M, K))
    P <- rep(0, K)
    for(j in 1:K) {
      sigma_inv[, , j] <- solve(sigma[, , j])
      P[j] <- mclust::dmvnorm(t(as.matrix(obs)),
                              mean = mu[, j],
                              sigma = sigma[, , j])
    }
    A <- (1:M)[index]
    B <- (1:M)[!index]
    # Yi Zhang, Gaussian Mixture Model Clustering with Incomplete Data (2021)
    xx1 <- fill_data_help1(obs = obs, B = B, mu = mu, alpha = p,
                           sigma_inv = sigma_inv, P = P)
    xx2 <- fill_data_help2(obs = obs, A = A, B = B, mu = mu, alpha = p,
                           sigma_inv = sigma_inv, P = P)
    obs[B] <- as.vector(xx1 %*% xx2)
  } 
  return(obs)
}

#' calculate the first half of the imputed values
fill_data_help1 <- function(obs, B, mu, alpha, sigma_inv, P) {
  K <- ncol(mu)
  l <- length(B)
  res <- matrix(rep(0, l * l), nrow = l)
  for(j in 1:K) {
    res <- res + alpha[j] * P[j] * sigma_inv[, , j][B, B]
  }
  return(solve(res))
}


#' calcualte the second half of the imputed values
fill_data_help2 <- function(obs, A, B, mu, alpha, sigma_inv, P) {
  K <- ncol(mu)
  l <- length(B)
  res <- rep(0, l)
  for(j in 1:K) {
    s1 <- sigma_inv[, , j][B, A, drop = FALSE]
    s2 <- sigma_inv[, , j][B, B, drop = FALSE]
    mu1 <- mu[A, j, drop = FALSE]
    mu2 <- mu[B, j, drop = FALSE]
    xx <- s1 %*% mu1 + s2 %*% mu2 - s1 %*% obs[A]
    res <- res + alpha[j] * P[j] * xx
  }
  return(res)
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
  cat("Variance-Covariance structure for biomarkers:", x$modelName, "model")
}