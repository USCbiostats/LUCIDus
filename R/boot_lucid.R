#' @title Inference of LUCID model based on bootstrap resampling
#' 
#' @description Generate \code{R} bootstrap replicates of LUCID parameters and 
#' derive confidence interval (CI) base on bootstrap. Bootstrap replicates are 
#' generated based on nonparameteric resampling, implemented by \code{ordinary} 
#' method of code{boot::boot} function.
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
#' @param model A LUCID model fitted by \code{est.lucid}.
#' @param conf A numeric scalar between 0 and 1 to specify confidence level(s) 
#' of the required interval(s).
#' @param R An integer to specify number of bootstrap replicates for LUCID model.
#' If feasible, it is recommended to set R > 1000. However, the convergence speed 
#' of LUCID varies greatly depending on data. If it takes very long time to run
#' 1000 replicates, it is recommend to set smaller values for R, such as 200.
#' 
#' @return A list, containing the following components:
#' \item{beta}{effect estimate for each exposure}
#' \item{mu}{cluster-specific mean for each omics feature}
#' \item{gamma}{effect estiamte for the association btween latent cluster and 
#' outcome}
#' \item{bootstrap}{The \code{boot} object returned by \code{boot:boot}}
#' 
#' @export
#' 
#' @import boot
#' @import progress
#' 
#' @examples
#' \dontrun{
#' # use simulated data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#' 
#' # fit lucid model
#' fit1 <- est_lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, 
#' seed = 1008)
#' 
#' # conduct bootstrap resampling
#' boot1 <- boot_lucid(G = G, Z = Z, Y = Y_normal, model = fit1, R = 100)
#' 
#' # check distribution for bootstrap replicates of the variable of interest
#' plot(boot1$bootstrap, 1)
#' 
#' # use 90% CI
#' boot2 <- boot_lucid(G = G, Z = Z, Y = Y_normal, model = fit1, R = 100, conf = 0.9)
#' }
boot_lucid <- function(G, 
                       Z, 
                       Y, 
                       CoG = NULL, 
                       CoY = NULL, 
                       model, 
                       conf = 0.95,
                       R = 100) {
  # prepare data for bootstrap (boot function require data in a matrix form, 
  # list data structure doesn't work)
  if(!is.null(model$selectG) | !is.null(model$selectZ)) {
    stop("Refit LUCID model with selected feature first then conduct bootstrap inference")
  }
  G <- as.matrix(G)
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)
  dimG <- ncol(G) 
  dimZ <- ncol(Z) 
  dimCoG <- ncol(CoG)
  dimCoY <- ncol(CoY)
  K <- model$K
  alldata <- as.data.frame(cbind(G, Z, Y, CoG, CoY))
  
  # bootstrap
  cat(paste0("Use Bootstrap resampling to derive ", 100 * conf, "% CI for LUCID \n"))
  #initialize progress bar object
  pb <- progress::progress_bar$new(total = R + 1) 
  bootstrap <- boot(data = alldata, 
                    statistic = lucid_par, 
                    R = R, 
                    dimG = dimG, 
                    dimZ = dimZ, 
                    dimCoY = dimCoY, 
                    dimCoG = dimCoG, 
                    model = model,
                    prog = pb)
  
  # bootstrap CIs
  ci <- gen_ci(bootstrap,
               conf = conf)
  
  # organize CIs
  beta <- ci[1:((K - 1) * dimG), ]
  mu <- ci[((K - 1) * dimG + 1): ((K - 1) * dimG + K * dimZ), ]
  gamma <- ci[-(1:((K - 1) * dimG + K * dimZ)), ]
  return(list(beta = beta,
              mu = mu, 
              gamma = gamma, 
              bootstrap = bootstrap))
}



# function to calculate parameters of LUCID model. use as statisitc input for
# boot function.
lucid_par <- function(data, indices, model, dimG, dimZ, dimCoY, dimCoG, prog) {
  #display progress with each run of the function
  prog$tick()
  Sys.sleep(0.01)
  
  # prepare data
  d <- data[indices, ]
  G <- as.matrix(d[, 1:dimG])
  Z <- as.matrix(d[, (dimG + 1):(dimG + dimZ)])
  Y <- as.matrix(d[, (dimG + dimZ + 1)])
  CoG <- CoY <- NULL
  K <- model$K
  if(!is.null(dimCoG)){
    CoG <- as.matrix(d[, (dimG + dimZ + 2):(dimG + dimZ + dimCoG + 1)])
  } 
  if(!is.null(dimCoY) && !is.null(dimCoG)){
    CoY <- as.matrix(d[, (dimG + dimZ + dimCoG + 1):ncol(d)])
  }
  if(!is.null(dimCoY) && is.null(dimCoG)){
    CoY <- as.matrix(d[, (dimG + dimZ + 2):ncol(d)])
  } 
  
  # fit lucid model
  seed <- sample(1:2000, 1)
  invisible(capture.output(try_lucid <- try(est_lucid(G = G, 
                                                      Z = Z, 
                                                      Y = Y,
                                                      CoY = CoY, 
                                                      CoG = CoG,
                                                      family = model$family, 
                                                      # modelName = model$modelName, 
                                                      modelName = NULL,
                                                      K = K, 
                                                      init_impute = model$init_impute,
                                                      init_par = model$init_par,
                                                      seed = seed))))
  if("try-error" %in% class(try_lucid)){
    n_par <- (K - 1) * dimG + K * dimZ + K
    if(!is.null(dimCoG)){
      n_par <- n_par + (K - 1) * dimCoG
    } 
    if(!is.null(dimCoY)){
      n_par <- n_par + dimCoY
    }
    par_lucid <- rep(0, n_par)
  } else{
    par_lucid <- c(as.vector(t(try_lucid$pars$beta)[-1, -1]),
                   as.vector(t(try_lucid$pars$mu)),
                   try_lucid$pars$gamma$beta)
    G_names <- as.vector(sapply(2:K, function(x) {
      paste0(colnames(try_lucid$pars$beta)[-1],
             ".cluster", x)
    }))
    Z_names <- as.vector(sapply(1:K, function(x) {
      paste0(colnames(try_lucid$pars$mu),
             ".cluster", x)
    }))
    if(is.null(names(try_lucid$pars$beta))) {
      Y_names <- paste0("cluster", 1:K)
    } else {
      Y_names <- names(try_lucid$pars$beta)
    }
    names(par_lucid) <- c(G_names, Z_names, Y_names)
    converge <- TRUE
  }
  return(par_lucid)
}


#' @title generate bootstrp ci (normal, basic and percentile)
#'
#' @param x an object return by boot function
#' @param conf A numeric scalar between 0 and 1 to specify confidence level(s) 
#' of the required interval(s).
#'
#' @return a matrix, the first column is t0 statistic from original model
#'
gen_ci <- function(x, conf = 0.95) {
  t0 <- x$t0
  res_ci <- NULL
  for (i in 1:length(t0)) {
    ci <- boot.ci(x, 
                  index = i, 
                  conf = conf, 
                  type = c("norm", "perc"))
    temp_ci <- c(ci$normal[2:3],
                 ci$percent[4:5])
    res_ci <- rbind(res_ci,
                    temp_ci)
  }
  res <- cbind(t0, res_ci)
  colnames(res) <- c("t0",
                     "norm_lower", "norm_upper",
                     "perc_lower", "perc_upper")
  return(res)
}
