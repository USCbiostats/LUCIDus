#' Fit a lucid model for integrated analysis on exposure, outcome and multi-omics data
#'
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
#' @param family Distribution of outcome. For continuous outcome, use "normal"; 
#' for binary outcome, use "binary". Default is "normal".
#' @param K Number of latent clusters (should be greater or equal than 2). 
#' Either an integer or a vector of integer. If K is a vector, model selection
#' on K is performed. 
#' @param Rho_G A scalar or a vector. This parameter is the LASSO penalty to regularize
#' exposures. If it is a vector, \code{lucid} will call \code{tune_lucid} to conduct
#' model selection and variable selection. User can try penalties from 0 to 1.
#' @param Rho_Z_Mu A scalar or a vector. This parameter is the LASSO penalty to 
#' regularize cluster-specific means for omics data (Z). If it is a vector, 
#' \code{lucid} will call \code{tune_lucid} to conduct model selection and 
#' variable selection. User can try penalties from 1 to 100.
#' @param Rho_Z_Cov A scalar or a vector. This parameter is the graphical LASSO
#' penalty to estimate sparse cluster-specific variance-covariance matrices for omics 
#' data (Z). If it is a vector, \code{lucid} will call \code{tune_lucid} to conduct
#' model selection and variable selection. User can try penalties from 0 to 1.
#' @param ... Other parameters passed to \code{est_lucid}
#'
#' @return An optimal lucid model
#' @export
#'
#' @examples
#' \dontrun{
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#' Y_binary <- sim_data$Y_binary
#' cov <- sim_data$Covariate
#' 
#' # fit lucid model
#' fit1 <- lucid(G = G, Z = Z, Y = Y_normal, family = "normal")
#' fit2 <- lucid(G = G, Z = Z, Y = Y_binary, family = "binary", useY = FALSE)
#' 
#' # including covariates
#' fit3 <- lucid(G = G, Z = Z, Y = Y_binary, family = "binary", CoG = cov)
#' fit4 <- lucid(G = G, Z = Z, Y = Y_binary, family = "binary", CoY = cov)
#' 
#' # tune K
#' fit5 <- lucid(G = G, Z = Z, Y = Y_binary, family = "binary", K = 2:5)
#' 
#' # variable selection
#' fit6 <- lucid(G = G, Z = Z, Y = Y_binary, family = "binary", Rho_G = seq(0.01, 0.1, by = 0.01))
#' fit7 <- lucid(G = G, Z = Z, Y = Y_binary, family = "binary", 
#' Rho_Z_Mu = seq(10, 100, by = 10), Rho_Z_Cov = 0.5,
#' init_par = "random", verbose_tune = TRUE)
#' }





lucid <- function(G, 
                  Z, 
                  Y, 
                  CoG = NULL, 
                  CoY = NULL, 
                  family = "normal", 
                  K = 2,
                  Rho_G = 0, 
                  Rho_Z_Mu = 0,
                  Rho_Z_Cov = 0, 
                  verbose_tune = FALSE,
                  ...) {
  # check data format
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
  
  # check K
  K <- as.integer(K)
  if(min(K) < 2) {
    stop("K should be integers greater than or equal to 2")
  }
  
  # check penalty terms
  check_Rho <- c(Rho_G, Rho_Z_Mu, Rho_Z_Cov)
  if(is.numeric(check_Rho)) {
    if(min(check_Rho) < 0) {
      stop("All penalties (Rho_G, Rho_Z_Mu and Rho_Z_Cov) should be greater than or equal to 0")
    }
  }
  
  flag_select_G <- FALSE
  flag_select_Z <- FALSE
  if(max(Rho_G) > 0) {
    flag_select_G <- TRUE
  }
  if(max(c(Rho_Z_Mu, Rho_Z_Cov)) > 0) {
    flag_select_Z <- TRUE
  }
  
  
  # tune lucid model
  if(verbose_tune) {
    res_tune <- tune_lucid(G = G, Z = Z, Y = Y, CoG = CoG, CoY = CoY,
                           family = family, K = K, 
                           Rho_G = Rho_G, Rho_Z_Mu = Rho_Z_Mu, Rho_Z_Cov = Rho_Z_Cov,
                           ...)
  } else {
    cat("Fitting LUCID model \n \n")
    invisible(capture.output(res_tune <- tune_lucid(G = G, Z = Z, Y = Y, CoG = CoG, CoY = CoY,
                                                    family = family, K = K, 
                                                    Rho_G = Rho_G, Rho_Z_Mu = Rho_Z_Mu, Rho_Z_Cov = Rho_Z_Cov,
                                                    ...)))
  }
  
  best_model <- res_tune$best_model
  
  # if variable selection is desirable, re-fit lucid model with selected variables
  select_G <- best_model$select$selectG
  if(flag_select_G) {
    if(sum(select_G) == 0) {
      cat("No exposure variables is selected using the given penalty Rho_G, please try a smaller one \n \n")
      cat("LUCID model will be fitted without variable selection on exposures (G) \n \n")
      select_G <- rep(TRUE, length(select_G))
    } else {
      cat(paste0(sum(select_G), "/", length(select_G)), "exposures are selected \n \n")
    }
  }
  
  select_Z <- best_model$select$selectZ
  if(flag_select_Z) {
    if(sum(select_Z) == 0) {
      cat("No omics variables is selected using the given penalty Rho_Z_Mu and Rho_Z_Cov, please try smaller ones \n \n")
      cat("LUCID model will be fitted without variable selection on omics data (Z) \n \n")
      select_Z <- rep(TRUE, length(select_Z))
    } else {
      cat(paste0(sum(select_Z), "/", length(select_Z)), "omics variables are selected \n \n")
      
    }
  }
  
  if(flag_select_G | flag_select_Z) {
    invisible(capture.output(best_model <- est_lucid(G = G[, select_G], Z = Z[, select_Z], Y = Y,
                                                     CoG = CoG, CoY = CoY, family = family,
                                                     K = K, ...)))
  }
  return(best_model) 
}