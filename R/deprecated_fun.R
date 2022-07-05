# This script documents deprecated functions. All the functions in this ducument
# should be deleted after release of next major update
lucid <- function(G, 
                       Z, 
                       Y, 
                       CoG = NULL, 
                       CoY = NULL, 
                       family = "normal", 
                       useY = TRUE,
                       K = 2:5,
                       Rho_G = 0, 
                       Rho_Z_Mu = 0,
                       Rho_Z_Cov = 0, 
                       ...) {
 stop("function `lucid` deprecated. Please use `tune_lucid`") 
}

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
 stop("function `est.lucid` deprecated. Please use `est_lucid`") 
}

boot.lucid <- function(G, 
                       Z, 
                       Y, 
                       CoG = NULL, 
                       CoY = NULL, 
                       model, 
                       conf = 0.95,
                       R = 100) {
  stop("function `boot.lucid` deprecated. Please use `boot_lucid`") 
}


