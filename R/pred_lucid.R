#' Predict the outcome based on a fitted LUCID model
#'
#' @param object A model fitted and returned by \code{\link{est.lucid}}
#' @param newG A new data set of genetic/environmental factors
#' @param newZ A new data set of biomarkers
#' @param newCoG Optional. A new data set of covariates included in the G->X analysis
#' @param newCoY Optional. A new data set of covariates included in the X->Y analysis
#' @param response Report the posterior distribution of cluster assignment (and the probability of binary outcome), default is TRUE
#' @param ... Other parameters to be passed to \code{predict}
#' @return A list contains predicted latent cluster and outcome for each observation
#' @export
#'
#' @examples
#' \dontrun{
#' index <- sample(1:2000, 200)
#' fit <- est.lucid(G = sim1[-index, 1:10], Z = sim1[-index, 11:20], Y = as.matrix(sim1[-index, 21]))
#' pred <- predict(object = fit, newG = sim1[index, 1:10], newZ = sim1[index, 11:20])
#' }

predict.lucid <- function(object, newG, newZ, newCoG = NULL, newCoY = NULL, response = TRUE, ...){
  n <- nrow(newG)
  gamma <- object$pars$gamma
  G <- cbind(newG, newCoG)
  Z <- newZ
  CoY <- newCoY
  ind.na <- Ind.NA(Z)
  K <- object$K
  dimCoY <- 0
  if(!is.null(CoY)){
    dimCoY <- ncol(CoY)
  }
  pars <- object$pars
  sigma.array <- array(as.numeric(unlist(pars$sigma)), dim = c(rep(ncol(newZ), 2), K))
  family.list <- switch(object$family, normal = normal(K = K, dimCoY), 
                        binary = binary(K = K, dimCoY))
  res <- Estep(beta = pars$beta, mu = pars$mu, sigma = sigma.array, gamma = pars$gamma,
               G = G, Z = Z, CoY = CoY, family.list, K = object$K, N = n, itr = 2, dimCoY = dimCoY, useY = FALSE, ind.na = ind.na)
  post.p <- res / rowSums(res) # posterior probability
  pred.x <- sapply(1:n, function(x) return(nnet::which.is.max(post.p[x, ])))
  e.cov <- rep(0, n)
  if(!is.null(newCoY)){
    e.cov <- newCoY %*% pars$gamma$beta[(K + 1):(K + ncol(newCoY))]
  }
  mu <- sapply(1:n, function(x) return(pars$gamma$beta[pred.x[x]] + e.cov[x]))
  if(object$family == "normal"){
    pred.y <- mu
  }
  if(object$family == "binary"){
    pred.y <- exp(mu) / (1 + exp(mu))
    if(response == TRUE){
      pred.y <- as.numeric(pred.y > 0.5)
    }
  }
  if(response == FALSE){
    pred.x <-  post.p
  }
  return(list(pred.x = pred.x, pred.y = pred.y))
}







