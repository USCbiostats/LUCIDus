# simulate multi-omics data for lucid
# simulation setting:
#

library(nnet)
library(mvtnorm)
library(mclust)


rm(list = ls())
source("./data-raw/sim_data_helper.R")


N <- 5000
K1 <- 2
K2 <- 3

nG <- 5
nZ1 <- 5
nZ2 <- 5

# step 1: generate G ====
set.seed(1008)
G <- t(replicate(N, rnorm(nG)))


# step 2: simulate X ====
# coef: G->X1
coef_GtoX1 <- matrix(rep(0, nG * K1), ncol = K1)
# leave the reference cluster as 0 (the 1st column) for convenience
coef_GtoX1[1:5, 2] <- c(log(3), log(2), log(1), log(1), log(1))
set.seed(1008)
X1 <- sim_X(G, coef_GtoX1)

# coef: G->X2
coef_GtoX2 <- matrix(rep(0, nG * K2), ncol = K2)
# leave the reference cluster as 0 (the 1st column) for convenience
coef_GtoX2[1:5, 2] <- c(log(3), log(2), log(1), log(1), log(1))
coef_GtoX2[1:5, 3] <- c(log(2), log(3), log(1), log(1), log(1))
set.seed(1008)
X2 <- sim_X(G, coef_GtoX2)


# step 3: simulate Z ====
# coef: X1->Z1
coef_X1toZ1 <- matrix(rep(0, nZ1 * K1), ncol = K1)
# edit the coef matrix, the column values represent the means for GMM.
coef_X1toZ1[1:5, 1] <- c(-3, -2, -1, 0, 0)
coef_X1toZ1[1:5, 2] <- c(3, 2, 1, 0, 0)
# edit the var-cov structure of Gaussian mixture model
sigma_X1toZ1 <- vector(mode = "list", length = K1)
sigma_X1toZ1[[1]] <- diag(nZ1)
sigma_X1toZ1[[2]] <- diag(nZ1)
set.seed(123)
Z1 <- sim_Z(X1, mu = coef_X1toZ1, sigma = sigma_X1toZ1)


# coef: X2->Z2
coef_X2toZ2 <- matrix(rep(0, nZ2 * K2), ncol = K2)
# edit the coef matrix, the column values represent the means for GMM.
coef_X2toZ2[1:5, 1] <- c(-3, -2, -1, 0, 0)
coef_X2toZ2[1:5, 2] <- c(-1, 0, 1, 0, 0)
coef_X2toZ2[1:5, 3] <- c(2, 1, 0, 0, 0)
# edit the var-cov structure of Gaussian mixture model
sigma_X2toZ2 <- vector(mode = "list", length = K2)
sigma_X2toZ2[[1]] <- diag(nZ2)
sigma_X2toZ2[[2]] <- diag(nZ2)
sigma_X2toZ2[[3]] <- diag(nZ2)
set.seed(123)
Z2 <- sim_Z(X2, mu = coef_X2toZ2, sigma = sigma_X2toZ2)



# step 4: simulate Y ====
set.seed(1008)
Cov1 <- rbinom(N, 1, 0.5) # sex
Cov2 <- rnorm(N, 1, 1) # 1st PC
Cov3 <- runif(N) # normalized age


# simulate a continuous Y without covariates
X <- data.frame(X1 = as.factor(X1),
                X2 = as.factor(X2))
b <- c(0, 1, -1, 2)
x <- model.matrix(~ X1 + X2, data = X)
mu_Y <- x %*% b
Y <- rnorm(N, mean = mu_Y)


## use regression to check results
fit1 <- multinom(X1 ~ G)
coef(fit1)
fit2 <- multinom(X2 ~ G)
coef(fit2)
fit3 <- Mclust(Z1)
fit3$parameters$mean
fit4 <- Mclust(Z2)
fit4$parameters$mean
fit5 <- lm(Y ~ X1 + X2, data = data.frame(Y = Y,
                                          X1 = as.factor(X1),
                                          X2 = as.factor(X2)))
coef(fit5)

sim_data_multi_omics <- list(G = G,
                             Z1 = Z1, 
                             Z2 = Z2,
                             Y = Y,
                             X1 = X1,
                             X2 = X2)
usethis::use_data(sim_data_multi_omics, overwrite = TRUE)
