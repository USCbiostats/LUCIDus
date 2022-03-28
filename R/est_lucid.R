#' @title  Fit LUCID model to conduct integrated clustering 
#' 
#' @description The Latent Unknown Clustering with Integrated Data (LUCID) performs 
#' integrative clustering using multi-view data. LUCID model is estimated via EM 
#' algorithm for model-based clustering. It also features variable selection,
#' integrated imputation, bootstrap inference and visualization via Sankey diagram.
#'
#' @param G Exposures, a numeric vector, matrix, or data frame. Categorical variable 
#' should be transformed into dummy variables. If a matrix or data frame, rows 
#' represent observations and columns correspond to variables.
#' @param Z Omics data, a numeric matrix or data frame. Rows correspond to observations
#' and columns correspond to variables.
#' @param Y Outcome, a numeric vector. Categorical variable is not allowed. Binary 
#' outcome should be coded as 0 and 1.
#' @param CoG Optional, covariates to be adjusted for estimating the latent cluster.
#' A numeric vector, matrix or data frame. Categorical variable should be transformed 
#' into dummy variables. 
#' @param CoY Optional, covariates to be adjusted for estimating the association 
#' between latent cluster and the outcome. A numeric vector, matrix or data frame. 
#' Categorical variable should be transformed into dummy variables.
#' @param K Number of latent clusters. An integer greater or equal to 2. User 
#' can use \code{\link{lucid}} to determine the optimal number of latent clusters.
#' @param family Distribution of outcome. For continuous outcome, use "normal"; 
#' for binary outcome, use "binary". Default is "normal".
#' @param useY Flag to include information of outcome when estimating the latent 
#' cluster. Default is TRUE.
#' @param tol Tolerance for convergence of EM algorithm. Default is 1e-3.
#' @param max_itr Max number of iterations for EM algorithm.
#' @param max_tot.itr Max number of total iterations for \code{est.lucid} function.
#' \code{est.lucid} may conduct EM algorithm for multiple times if the algorithm 
#' fails to converge.
#' @param Rho_G A scalar. Penalty to conduct LASSO regularization and obtain a sparse estimation
#' for effect of exposures. If user wants to tune the penalty, use the wrapper 
#' function \code{lucid}
#' @param Rho_Z_Mu A scalar. Penalty to conduct LASSO regularization and obtain a sparse 
#' estimation of cluster-specific mean for omics data. If user wants to tune the 
#' penalty, use the wrapper function \code{lucid}
#' @param Rho_Z_Cov A scalar. Penalty to conduct graphic LASSO regularization and obtain a
#' sparse estimation of cluster-specific variance-covariance matrices for omics 
#' data. If user wants to tune the penalty, use the wrapper function \code{lucid}
#' @param modelName The variance-covariance structure for omics data. 
#' See \code{mclust::mclustModelNames} for details.
#' @param seed An integer to initialize the EM algorithm or imputing missing values.
#'Default is 123.
#' @param init_impute Method to initialize the imputation of missing values in 
#' LUCID. "mclust" will use \code{mclust:imputeData} to implement EM Algorithm 
#' for Unrestricted General Location Model to impute the missing values in omics 
#' data; \code{lod} will initialize the imputation via relacing missing values by
#' LOD / sqrt(2). LOD is determined by the minimum of each variable in omics data.
#' @param init_par Method to initialize the EM algorithm. "mclust" will use mclust
#' model to initialize parameters; "random" initialize parameters from uniform 
#' distribution.
#' @param verbose A flag indicates whether detailed information for each iteration
#' of EM algorithm is printed in console. Default is FALSE.
#' 
#' 
#' 
#' @return A list which contains the several features of LUCID, including:
#' \item{pars}{Estimates of parameters of LUCID, including beta (effect of 
#' exposure), mu (cluster-specific mean for omics data), sigma (cluster-specific 
#' variance-covariance matrix for omics data) and gamma (effect estimate of association
#' between latent cluster and outcome)}
#' \item{K}{Number of latent cluster}
#' \item{modelName}{Geometric model to estiamte variance-covariance matrix for 
#' omics data}
#' \item{likelihood}{The log likelihood of the LUCID model}
#' \item{post.p}{Posterior inclusion probability (PIP) for assigning observation i
#' to latent cluster j}
#' \item{Z}{If missing values are observed, this is the complet dataset for omics
#' data with missing values imputed by LUCID}
#' 
#' @importFrom nnet multinom
#' @import mclust
#' @importFrom glmnet glmnet
#' @importFrom glasso glasso
#' @import stats
#' @import utils
#' @export
#'
#' @references
#' Cheng Peng, Jun Wang, Isaac Asante, Stan Louie, Ran Jin, Lida Chatzi, 
#' Graham Casey, Duncan C Thomas, David V Conti, A Latent Unknown Clustering 
#' Integrating Multi-Omics Data (LUCID) with Phenotypic Traits, Bioinformatics,
#' btz667, https://doi.org/10.1093/bioinformatics/btz667.
#' 
#' 
#' @examples
#' \dontrun{
#' # use simulated data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#' Y_binary <- sim_data$Y_binary
#' cov <- sim_data$Covariate
#' 
#' # fit LUCID model with continuous outcome
#' fit1 <- est.lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, 
#' seed = 1008)
#' 
#' # fit LUCID model with block-wise missing pattern in omics data
#' Z_miss_1 <- Z
#' Z_miss_1[sample(1:nrow(Z), 0.3 * nrow(Z)), ] <- NA
#' fit2 <- est.lucid(G = G, Z = Z_miss_1, Y = Y_normal, family = "normal", K = 2)
#' 
#' # fit LUCID model with sporadic missing pattern in omics data
#' Z_miss_2 <- Z
#' index <- arrayInd(sample(length(Z_miss_2), 0.3 * length(Z_miss_2)), dim(Z_miss_2))
#' Z_miss_2[index] <- NA
#' # initialize imputation by imputing 
#' fit3 <- est.lucid(G = G, Z = Z_miss_2, Y = Y_normal, family = "normal", 
#' K = 2, seed = 1008, init_impute = "lod") 
#' LOD
#' # initialize imputation by mclust
#' fit4 <- est.lucid(G = G, Z = Z_miss_2, Y = Y, family = "normal", K = 2, 
#' seed = 123, init_impute = "mclust") 
#' 
#' # fit LUCID model with binary outcome
#' fit5 <- est.lucid(G = G, Z = Z, Y = Y_binary, family = "binary", K = 2,
#' seed = 1008)
#' 
#' # fit LUCID model with covariates
#' fit6 <- est.lucid(G = G, Z = Z, Y = Y_binary, CoY = cov, family = "binary", 
#' K = 2, seed = 1008)
#' 
#' # use LUCID model to conduct integrated variable selection
#' # select exposure
#' fit6 <- est.lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal", 
#' K = 2, seed = 1008, Rho_G = 0.1)
#' # select omics data
#' fit7 <- est.lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal",
#' K = 2, seed = 1008, Rho_Z_Mu = 90, Rho_Z_Cov = 0.1, init_par = "random")
#' 
#' }
est.lucid <- function(G, 
                      Z, 
                      Y, 
                      CoG = NULL, 
                      CoY = NULL, 
                      K = 2, 
                      family = c("normal", "binary"), 
                      useY = TRUE, 
                      tol = 1e-3,
                      max_itr = 1e3,
                      max_tot.itr = 1e4,
                      Rho_G = 0,
                      Rho_Z_Mu = 0,
                      Rho_Z_Cov = 0, 
                      modelName = "VVV",
                      seed = 123,
                      init_impute = c("mclust", "lod"),
                      init_par = c("mclust", "random"),
                      verbose = FALSE) {
  
  # 1. basic setup for estimation function =============
  family <- match.arg(family)
  init_impute <- match.arg(init_impute)
  init_par <- match.arg(init_par)
  Select_G <- FALSE
  Select_Z <- FALSE
  if(Rho_G != 0) {
    Select_G <- TRUE
  }
  if(Rho_Z_Mu != 0 | Rho_Z_Cov != 0) {
    Select_Z <- TRUE
  }
  
  
  ## 1.1 check data format ====
  if(is.null(G)) {
    stop("Input data 'G' is missing")
  } else {
    if(!is.matrix(G)) {
      G <- as.matrix(G)
      if(!is.numeric(G)) {
        stop("Input data 'G' should be numeric; categorical variables should be transformed into dummies")
      }
    }
  }
  if(is.null(colnames(G))){
    Gnames <- paste0("G", 1:ncol(G))
  } else {
    Gnames <- colnames(G)
  }
  colnames(G) <- Gnames
  
  if(is.null(Z)) {
    stop("Input data 'Z' is missing")
  } else {
    if(!is.matrix(Z)) {
      Z <- as.matrix(Z)
      if(!is.numeric(Z)) {
        stop("Input data 'Z' should be numeric")
      }
    }
  }
  if(is.null(colnames(Z))){
    Znames <- paste0("Z", 1:ncol(Z))
  } else {
    Znames <- colnames(Z)
  }
  
  if(is.null(Y)) {
    stop("Input data 'Y' is missing")
  } else {
    if(!is.matrix(Y)) {
      Y <- as.matrix(Y)
      if(!is.numeric(Y)) {
        stop("Input data 'Y' should be numeric; binary outcome should be transformed them into dummies")
      }
      if(ncol(Y) > 1) {
        stop("Only continuous 'Y' or binary 'Y' is accepted")
      }
    }
  }
  if(is.null(colnames(Y))) {
    Ynames <- "outcome"
  } else {
    Ynames <- colnames(Y)
  }
  colnames(Y) <- Ynames
  if(family == "binary") {
    if(!(all(Y %in% c(0, 1)))) {
      stop("Binary outcome should be coded as 0 and 1")
    }
  }
  
  
  CoGnames <- NULL
  if(!is.null(CoG)) {
    if(!is.matrix(CoG)) {
      CoG <- as.matrix(CoG)
      if(!is.numeric(CoG)) {
        stop("Input data 'CoG' should be numeric; categroical variables should be transformed into dummies")
      }
    }
    if(is.null(colnames(CoG))) {
      CoGnames <- paste0("CoG", 1:ncol(CoG))
    } else {
      CoGnames <- colnames(CoG)  
    }
    colnames(CoG) <- CoGnames
  }
  
  CoYnames <- NULL
  if(!is.null(CoY)) {
    if(!is.matrix(CoY)) {
      CoY <- as.matrix(CoY)
      if(!is.numeric(CoY)) {
        stop("Input data 'CoY' should be numeric; categorical variables should be transformed into dummies")
      }
    }
    if(is.null(colnames(CoY))) {
      CoYnames <- paste0("CoY", 1:ncol(CoY))
    } else {
      CoYnames <- colnames(CoY)
    }
    colnames(CoY) <- CoYnames
  }
  
  ## 1.2 record input dimensions, family function ====
  N <- nrow(Y)
  dimG <- ncol(G)
  dimZ <- ncol(Z); 
  dimCoG <- ifelse(is.null(CoG), 0, ncol(CoG))
  dimCoY <- ifelse(is.null(CoY), 0, ncol(CoY))
  G <- cbind(G, CoG)
  Gnames <- c(Gnames, CoGnames)
  family.list <- switch(family, normal = normal(K = K, dimCoY), 
                        binary = binary(K = K, dimCoY))
  Mstep_Y <- family.list$f.maxY
  switch_Y <- family.list$f.switch
  
  
  ## 1.3. check missing pattern ====
  na_pattern <- check_na(Z)
  if(na_pattern$impute_flag) {
    # initialize imputation
    if(init_impute == "mclust") {
      cat("Intializing imputation of missing values in 'Z' via the mix package \n\n")
      invisible(capture.output(Z <- mclust::imputeData(Z, seed = seed)))
      Z[na_pattern$indicator_na == 3, ] <- NA  
    }
    if(init_impute == "lod") {
      cat("Intializing imputation of missing values in 'Z' via LOD / sqrt(2) \n\n")
      Z <- apply(Z, 2, fill_data_lod)
      colnames(Z) <- Znames
    }
    
  }
  
  
  # 2. EM algorithm for LUCID ================
  tot.itr <- 0
  convergence <- FALSE
  while(!convergence && tot.itr <= max_tot.itr) {
    if(tot.itr > 0) {
      seed <- seed + 10
    }
    set.seed(seed)
    
    ## 2.1 initialize model parameters ====
    
    # initialize beta 
    res.beta <- matrix(data = runif(K * (dimG + dimCoG + 1)), nrow = K) 
    res.beta[1, ] <- 0
    
    # initialize mu and sigma
    # initialize by mclust
    if(init_par == "mclust") {
      cat("Initialize LUCID with mclust \n")
      invisible(capture.output(mclust.fit <- Mclust(Z[na_pattern$indicator_na != 3, ], 
                                                    G = K,
                                                    modelNames = modelName)))
      if(is.null(mclust.fit)) {
        stop("mclust failed for specified model - please set modelNames to `NULL` to conduct automatic model selection ")
      }
      if(is.null(modelName)){
        model.best <- mclust.fit$modelName
      } else{
        model.best <- modelName
      }
      res.mu <- t(mclust.fit$parameters$mean)
      res.sigma <- mclust.fit$parameters$variance$sigma
    } else { # initialize by random guess
      cat("Initialize LUCID with random values from uniform distribution \n")
      if(is.null(modelName)){
        model.best <- "VVV"
        warning("GMM model for LUCID is not specified, 'VVV' model is used by default")
      } else{
        model.best <- modelName
      }
      res.mu <- matrix(runif(dimZ * K, min = -0.5, max = 0.5),
                       nrow = K)
      
      res.sigma <- gen_cov_matrices(dimZ = dimZ, K = K)
    }
    
    # initialize family specific parameters gamma
    res.gamma <- family.list$initial.gamma(K, dimCoY)
    
    
    # start EM algorithm 
    cat(paste0("Fitting LUCID model (K = ", K, ") \n"))
    res.loglik <- -Inf
    itr <- 0
    while(!convergence && itr <= max_itr){
      itr <- itr + 1
      tot.itr <- tot.itr + 1
      check.gamma <-  TRUE
      
      # 2.2 E-step ====
      # calculate log-likelihood for observation i being assigned to cluster j
      new.likelihood <- Estep(beta = res.beta, 
                              mu = res.mu, 
                              sigma = res.sigma, 
                              gamma = res.gamma,
                              G = G, 
                              Z = Z, 
                              Y = Y, 
                              CoY = CoY, 
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
        cat("iteration", itr,": EM algorithm collapsed: invalid estiamtes due to over/underflow, try LUCID with another seed \n")
        break
      } else{
        if(isTRUE(verbose)) {
          cat("iteration", itr,": E-step finished.\n")  
        }
      }
      
      
      # 2.3 M-step - parameters ====
      # update model parameters to maximize the expected likelihood
      invisible(capture.output(new.beta <- Mstep_G(G = G, 
                                                   r = res.r, 
                                                   selectG = Select_G, 
                                                   penalty = Rho_G, 
                                                   dimG = dimG, 
                                                   dimCoG = dimCoG,
                                                   K = K)))
      new.mu.sigma <- Mstep_Z(Z = Z, 
                              r = res.r, 
                              selectZ = Select_Z, 
                              penalty.mu = Rho_Z_Mu, 
                              penalty.cov = Rho_Z_Cov,
                              model.name = model.best, 
                              K = K, 
                              ind.na = na_pattern$indicator_na, 
                              mu = res.mu)
      if(is.null(new.mu.sigma$mu)){
        cat("variable selection failed, try LUCID with another seed \n")
        break
      }
      if(useY){
        new.gamma <- Mstep_Y(Y = Y, r = res.r, CoY = CoY, K = K, CoYnames)
        check.gamma <- is.finite(unlist(new.gamma))
      }
      
      
      # 2.4 M step - impute missing values ====
      if(na_pattern$impute_flag){
        Z <- Istep_Z(Z = Z, 
                     p = res.r, 
                     mu = res.mu, 
                     sigma = res.sigma,
                     index = na_pattern$index)
      }
      
      
      # 2.5 control step ====
      check.value <- all(is.finite(new.beta), 
                         is.finite(unlist(new.mu.sigma)), 
                         check.gamma)
      
      if(!check.value){
        cat("iteration", itr,": Invalid estimates, try LUCID with another seed \n")
        break
      } else{
        res.beta <- new.beta
        res.mu <- new.mu.sigma$mu
        res.sigma <- new.mu.sigma$sigma
        if(useY){
          res.gamma <- new.gamma
        }

        new.loglik <- sum(rowSums(res.r * new.likelihood))

        if(Select_G) {
          new.loglik <- new.loglik - Rho_G * sum(abs(res.beta))
        }
        if(Select_Z) {
          new.loglik <- new.loglik - Rho_Z_Mu * sum(abs(res.mu)) - Rho_Z_Cov * sum(abs(res.sigma))
        }
        if(isTRUE(verbose)) {
          if(Select_G | Select_Z) {
            cat("iteration", itr,": M-step finished, ", "penalized loglike = ", sprintf("%.3f", new.loglik), "\n")
          } else{
            cat("iteration", itr,": M-step finished, ", "loglike = ", sprintf("%.3f", new.loglik), "\n")
          }
        } else {
          cat(".")
        }
        
        if(abs(res.loglik - new.loglik) < tol){
          convergence <- TRUE
          cat("Success: LUCID converges!", "\n\n")
        }
        res.loglik <- new.loglik
      }
    }
  }
  
  # 3. summarize results ===============
  if(!useY){
    res.gamma <- Mstep_Y(Y = Y, r = res.r, CoY = CoY, K = K, CoYnames = CoYnames)
  }
  # browser()
  res.likelihood <- Estep(beta = res.beta, 
                          mu = res.mu, 
                          sigma = res.sigma, 
                          gamma = res.gamma,
                          G = G, 
                          Z = Z, 
                          Y = Y, 
                          family.list = family.list, 
                          itr = itr, 
                          CoY = CoY, 
                          N = N, 
                          K = K, 
                          dimCoY = dimCoY, 
                          useY = useY, 
                          ind.na = na_pattern$indicator_na)
  res.r <- t(apply(res.likelihood, 1, lse_vec))
  

  res.loglik <- sum(rowSums(res.r * res.likelihood))

  if(Select_G) {
    res.loglik <- res.loglik - Rho_G * sum(abs(res.beta))
  }
  if(Select_Z) {
    res.loglik <- res.loglik - Rho_Z_Mu * sum(abs(res.mu)) - Rho_Z_Cov * sum(abs(res.sigma))
  }
  pars <- switch_Y(beta = res.beta, mu = res.mu, sigma = res.sigma, gamma = res.gamma, K = K)
  res.r <- res.r[, pars$index]
  colnames(pars$beta) <- c("intercept", Gnames)
  colnames(pars$mu) <- Znames
  if(Select_G){
    tt1 <- apply(pars$beta[, -1], 2, range)
    selectG <- abs(tt1[2, ] - tt1[1, ]) > 0.001
  } else{
    selectG <- rep(TRUE, dimG)
  }
  if(Select_Z){
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
                  select = list(selectG = selectG, selectZ = selectZ), 
                  useY = useY,
                  Z = Z,
                  init_impute = init_impute,
                  init_par = init_par)
  class(results) <- c("lucid")
  return(results)
}