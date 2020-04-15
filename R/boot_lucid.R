boot.lucid <- function(G, Z, Y, CoG = NULL, CoY = NULL, model, R = 100){
  dimG <- ncol(G); dimZ <- ncol(Z); dimCoY <- ncol(CoY); dimCoG  <- ncol(CoG); K <- model$K
  alldata <- as.data.frame(cbind(G, Z, Y, CoG, CoY))
  bootstrap <- boot(data = alldata, statistic = lucid_par, R = R, parallel = "multicore", ncpus = detectCores(),
                    dimG = dimG, dimZ = dimZ, dimCoY = dimCoY, dimCoG = dimCoG, model = model)
  sd <- sapply(1:length(bootstrap$t0), function(x) sd(bootstrap$t[, x]))
  model.par <- c(model$pars$beta[-1, -1], as.vector(t(model$pars$mu)), model$pars$gamma$beta)
  dd <- data.frame(original = model.par,
                   lower = model.par - 1.96 * sd, upper = model.par + 1.96 * sd)
  return(dd)
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
                               family = model$family, control = model$par.control, tune = model$par.tune,
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