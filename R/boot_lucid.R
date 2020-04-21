#' @title Bootstrap 
#' 
#' @description This function provides SEs of parameter estimates from a LUCID model through bootstrap method.
#'
#' @param G Genetic features/environmental exposures, a \code{\link{matrix}}.
#' @param Z Biomarkers/other omics data, a \code{\link{matrix}}.
#' @param Y Disease outcome, it is suggested to transform it into a n by 1 \code{\link{matrix}}.
#' @param CoG Optional, matrix. Covariates to be adjusted for estimating the latent cluster.
#' @param CoY Optional, matrix. Covariates to be adjusted for estimating the outcome.
#' @param model A LUCID model fitted by \code{\link{est.lucid}}.
#' @param R Number of boostrap iterations.
#' @return A list of estimates with thier 95 percent CI.
#' \item{Beta}{Estimates, SE and 95% CI of genetic effects/environmental exposure (and effects of covariates, if included), matrix}
#' \item{Mu}{Estimates, SE and 95% CI of cluster-specific biomarker means, matrix}
#' \item{Gamma}{Estimates, SE and 95% CI of cluster-specific disease risk (and effects of covariates, if included), matrix}
#' @export
#' @import boot
#' @import parallel
#' @author Yinqi Zhao, Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Cheng Peng, Jun Wang, Isaac Asante, Stan Louie, Ran Jin, Lida Chatzi, Graham Casey, Duncan C Thomas, David V Conti, A Latent Unknown Clustering Integrating Multi-Omics Data (LUCID) with Phenotypic Traits, Bioinformatics, , btz667, https://doi.org/10.1093/bioinformatics/btz667.
#' @examples
#' 
boot.lucid <- function(G, Z, Y, CoG = NULL, CoY = NULL, model, R = 100){
  ss <- model$select
  G <- G[, ss$selectG]
  Z <- Z[, ss$selectZ]
  dimG <- ncol(G); dimZ <- ncol(Z); dimCoY <- ncol(CoY); dimCoG  <- ncol(CoG); K <- model$K
  alldata <- as.data.frame(cbind(G, Z, Y, CoG, CoY))
  bootstrap <- boot(data = alldata, statistic = lucid_par, R = R, parallel = "multicore", ncpus = detectCores(),
                    dimG = dimG, dimZ = dimZ, dimCoY = dimCoY, dimCoG = dimCoG, model = model)
  sd <- sapply(1:length(bootstrap$t0), function(x) sd(bootstrap$t[, x]))
  model.par <- c(model$pars$beta[-1, c(FALSE, ss$selectG)], as.vector(t(model$pars$mu[, ss$selectZ])), model$pars$gamma$beta)
  dd <- data.frame(original = model.par, sd = sd, 
                   lower = model.par - 1.96 * sd, upper = model.par + 1.96 * sd)
  beta <- dd[1:((K - 1) * dimG), ]
  row.names(beta) <- paste0(colnames(G), ".cluster", 2:K)
  mu <- dd[((K - 1) * dimG + 1): ((K - 1) * dimG + K * dimZ), ]
  row.names(mu) <- paste0(rep(colnames(Z), K), ".cluster", sapply(1:K, function(x) rep(x, dimZ)))
  gamma <- dd[-(1:((K - 1) * dimG + K * dimZ)), ]
  row.names(gamma) <- c("reference", paste0("cluster", 2:K), colnames(CoY))
  return(structure(list(beta = beta, mu = mu, gamma = gamma)))
}

lucid_par <- function(data, indices, dimG, dimZ, dimCoY, dimCoG, model) {
  d <- data[indices, ]
  G <- as.matrix(d[, 1:dimG])
  Z <- as.matrix(d[, (dimG + 1):(dimG + dimZ)])
  Y <- as.matrix(d[, (dimG + dimZ + 1)])
  CoG <- CoY <- NULL
  if(!is.null(dimCoG)){
    CoG <- as.matrix(d[, (dimG + dimZ + 2):(dimG + dimZ + dimCoG + 1)])
  } 
  if(!is.null(dimCoY) && !is.null(dimCoG)){
    CoY <- as.matrix(d[, (dimG + dimZ + dimCoG + 1):ncol(d)])
  }
  if(!is.null(dimCoY) && is.null(dimCoG)){
    CoY <- as.matrix(d[, (dimG + dimZ + 2):ncol(d)])
  } 
  converge <- FALSE
  while(!converge){
    try_lucid <- try(est.lucid(G = G, 
                               Z = Z, 
                               Y = Y,
                               CoY = CoY, 
                               CoG = CoG,
                               family = model$family, control = model$par.control,
                               Z.var.str = model$Z.var.str, K = model$K))
    if("try-error" %in% class(try_lucid)){
      next
    } else{
      par_lucid <- c(try_lucid$pars$beta[-1, -1],
                     as.vector(t(try_lucid$pars$mu)),
                     try_lucid$pars$gamma$beta)
      converge <- TRUE
    }
  }
  return(par_lucid)
}