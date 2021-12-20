#' Check missing patterns in omics data Z
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



#' I-step of LUCID: impute missing data in Z by maximizing the likelihood
#'
#' @param Z 
#' @param r 
#' @param mu
#' @param sigma
#' @param index
#'
#' @return a complete dataset of Z
#' 
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

#' Calculate the first half of the imputed values
fill_data_help1 <- function(obs, B, mu, alpha, sigma_inv, P) {
  K <- ncol(mu)
  l <- length(B)
  res <- matrix(rep(0, l * l), nrow = l)
  for(j in 1:K) {
    res <- res + alpha[j] * P[j] * sigma_inv[, , j][B, B]
  }
  return(solve(res))
}


#' Calcualte the second half of the imputed values
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
