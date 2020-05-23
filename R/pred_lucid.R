#' Predict the outcome based on a fitted LUCID model
#'
#' @param object A model fitted and returned by \code{\link{est.lucid}}
#' @param newG A new data set of genetic/environmental factors
#' @param newZ A new data set of biomarkers
#' @param newCoG Optional. A new data set of covariates included in the G->X analysis
#' @param newCoY Optional. A new data set of covariates included in the X->Y analysis
#' @param ... Other parameters to be passed to \code{predict}
#' @return A list contains predicted latent cluster and outcome for each observation
#' @export
#'
#' @examples
#' \dontrun{
#' index <- sample(1:3000, 200)
#' fit <- est.lucid(G = G1[-index, ], Z = Z1[-index, ], Y = as.matrix(Y1[-index, ]))
#' pred <- predict(object = fit, newG = G1[index, ], newZ = Z1[index, ])
#' }

predict.lucid <- function(object, newG, newZ, newCoG = NULL, newCoY = NULL, ...){
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
  post.p <- res / rowSums(res) # posterior probablity
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
  }
  return(list(pred.x = pred.x, pred.y = pred.y))
}







