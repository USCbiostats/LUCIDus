# codes to prepare simulated dataset ====

source("./data-raw/sim_data_helper.R")
library(mvtnorm)
N <- 2000
K <- 2
nG <- 10
nZ <- 10

## step 1: generate G====
set.seed(1008)
G <- t(replicate(N, rnorm(nG)))

## step 2: simulate X ====
# coef: G->X
coef_GtoX <- matrix(rep(0, nG * K), ncol = K)
# edit the coef matrix
# leave the reference cluster as 0 (the 1st column) for convenience
coef_GtoX[1:4, 2] <- c(log(2), log(2), log(2), log(2))
X <- sim_X(G, coef_GtoX)


## step 3: simulate Z ====
# coef: X->Z
coef_XtoZ <- matrix(rep(0, nZ * K), ncol = K)
# edit the coef matrix, the column values represent the means for GMM.
coef_XtoZ[1:5, 1] <- rep(-1, 5)
coef_XtoZ[1:5, 2] <- rep(1, 5)
# edit the var-cov structure of Gaussian mixture model
sigma_XtoZ <- vector(mode = "list", length = K)
t1 <- matrix(runif(nZ^2, min = -0.5, max = 0.5), nrow = nZ)
sigma_XtoZ[[1]] <- t1 %*% t(t1)
t2 <- matrix(runif(nZ^2, min = -0.5, max = 0.5), nrow = nZ)
sigma_XtoZ[[2]] <- t2 %*% t(t2)
set.seed(123)
Z <- sim_Z(X, mu = coef_XtoZ, sigma = sigma_XtoZ)


## step 4: simulate Y ====

### 4.1: continuous Y ====
CovX <- NULL
# beta matrix: the first k elements represent the mean for each cluster
coef_XtoY <- c(0.5, 1.5)
set.seed(1008)
Y_normal <- as.matrix(sim_Y_normal(X = X, 
                                   CovX = CovX, 
                                   beta = coef_XtoY, 
                                   sigma = c(1, 1)))


# test the simulated data
# G <- as.matrix(G)
# Z <- as.matrix(Z)
# Y <- as.matrix(Y)
# fit1 <- est.lucid(G = G, 
#                   Z = Z, 
#                   Y = Y_normal, 
#                   family = "normal", 
#                   K = 2,
#                   seed = 1008)
# summary(fit1)

colnames(G) <- paste0("exposure", 1:10)
colnames(Z) <- paste0("metabolite", 1:10)
colnames(Y_normal) <- "pfas_concentration"


### 4.2: binary Y ====
CovX <- t(replicate(N, rnorm(2)))

# beta matrix: the first k elements represent the mean for each cluster
coef_XtoY <- c(0, log(2.2), -0.2, 0.2)
set.seed(1008)
Y_binary <- as.matrix(sim_Y_binary(X = X, CovX = CovX, beta = coef_XtoY))


# test the simulated data
# G <- as.matrix(G)
# Z <- as.matrix(Z)
# Y <- as.matrix(Y)
# fit2 <- est.lucid(G = G,
#                   Z = Z,
#                   Y = Y_binary,
#                   CoY = CovX,
#                   family = "binary",
#                   K = 2,
#                   seed = 1008)
# summary(fit2)

colnames(Y_binary) <- "liver_injury"
colnames(CovX) <- c("Cov1", "Cov2")

sim_data <- list(G = G,
                 Z = Z,
                 Y_normal = Y_normal,
                 Y_binary = Y_binary,
                 Covariate = CovX,
                 X = X)
usethis::use_data(sim_data, overwrite = TRUE)
