#' @title Predict cluster assignment and outcome based on LUCID model
#'
#' @param model A model fitted and returned by \code{\link{est.lucid}}
#' @param G A new data set of genetic/environmental factors
#' @param Z A new data set of biomarkers
#' @param CoG Optional. A new data set of covariates included in the G->X analysis
#' @param CoY Optional. A new data set of covariates included in the X->Y analysis
#' @param response Report the posterior distribution of cluster assignment (and the probability of binary outcome), default is TRUE
#' @param ... Other parameters to be passed to \code{predict}
#' @return A list contains predicted latent cluster and outcome for each observation
#' @export
#'
#' @examples
#' \dontrun{
#' # prepare data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#' 
#' # fit lucid model
#' fit1 <- est.lucid(G = G, Z = Z, Y = Y_normal, K = 2, family = "normal")
#' 
#' # prediction on training set
#' pred1 <- predict_lucid(model = fit1, G = G, Z = Z)
#' }

predict_lucid <- function(model, 
                          G, 
                          Z, 
                          CoG = NULL, 
                          CoY = NULL, 
                          response = TRUE, ...){
  
  if(class(model) != "lucid") {
    stop("model should be fitted by est.lucid")
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
  
  if(!is.null(CoG)) {
    if(!is.matrix(CoG)) {
      CoG <- as.matrix(CoG)
      if(!is.numeric(CoG)) {
        stop("Input data 'CoG' should be numeric; categroical variables should be transformed into dummies")
      }
    }
  }
  
  if(!is.null(CoY)) {
    if(!is.matrix(CoY)) {
      CoY <- as.matrix(CoY)
      if(!is.numeric(CoY)) {
        stop("Input data 'CoY' should be numeric; categorical variables should be transformed into dummies")
      }
    }
  }
  
  n <- nrow(G)
  K <- model$K
  
  # model parameters
  beta <- model$pars$beta
  mu <- model$pars$mu
  Sigma <- model$pars$sigma
  Sigma.array <- array(as.numeric(unlist(Sigma)), dim = c(rep(ncol(Z), 2), K))
  gamma <- model$pars$gamma
  
  G <- cbind(G, CoG)
  na_pattern <- check_na(Z)
  dimCoY <- 0
  if(!is.null(CoY)){
    dimCoY <- ncol(CoY)
  }
  family.list <- switch(model$family, 
                        normal = normal(K = K, dimCoY), 
                        binary = binary(K = K, dimCoY))
  
  res <- Estep(beta = beta, 
               mu = mu, 
               sigma = Sigma.array, 
               gamma = gamma,
               G = G, 
               Z = Z, 
               CoY = CoY, 
               family.list, 
               K = K, 
               N = n,
               itr = 2,
               dimCoY = dimCoY, 
               useY = FALSE, 
               ind.na = na_pattern$indicator_na)
  
  # normalize the log-likelihood to probability
  res.r <- t(apply(res, 1, lse_vec))
  # predicted latent cluster
  pred.x <- sapply(1:n, function(x) return(nnet::which.is.max(res.r[x, ])))
  mu_Y <- cbind(res.r, CoY) %*% gamma$beta
  
  
  if(model$family == "normal"){
    pred.y <- mu_Y
  }
  if(model$family == "binary"){
    pred.y <- exp(mu_Y) / (1 + exp(mu_Y))
    if(response == TRUE){
      pred.y <- as.numeric(pred.y > 0.5)
    }
  }
  
  return(list(post.p = res.r,
              pred.x = pred.x, 
              pred.y = pred.y))
}







