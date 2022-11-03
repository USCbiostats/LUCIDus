## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(LUCIDus)

## ----out.width="70%", echo=FALSE----------------------------------------------
knitr::include_graphics("DAG.png")

## ----out.width="70%", echo=FALSE----------------------------------------------
knitr::include_graphics("workflow.png")

## ---- eval=FALSE--------------------------------------------------------------
#  library(LUCIDus)
#  # use simulated data
#  G <- sim_data$G
#  Z <- sim_data$Z
#  Y_normal <- sim_data$Y_normal
#  Y_binary <- sim_data$Y_binary
#  cov <- sim_data$Covariate
#  
#  # fit LUCID model with continuous outcome
#  fit1 <- lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, seed = 1008)
#  
#  # fit LUCID model with binary outcome
#  fit2 <- lucid(G = G, Z = Z, Y = Y_binary, family = "binary", K = 2, seed = 1008)
#  
#  # fit LUCID model with covariates
#  fit3 <- lucid(G = G, Z = Z, Y = Y_binary, CoY = cov, family = "binary", K = 2, seed = 1008)
#  fit4 <- lucid(G = G, Z = Z, Y = Y_binary, CoG = cov, family = "binary", K = 2, seed = 1008)

## ---- eval=FALSE--------------------------------------------------------------
#  # unsupervised lucid
#  fit5 <- lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, useY = FALSE, seed = 1008)

## ---- eval=FALSE--------------------------------------------------------------
#  # fit LUCID model with automatic selection on covariance models
#  fit6 <- lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, modelName = NULL, seed = 1008)
#  # check the optimal model
#  fit6$modelName
#  
#  # fit LUCID model with a specified covariance model
#  fit7 <- lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, modelName = "EII", seed = 1008)

## ---- eval=FALSE--------------------------------------------------------------
#  # initialize EM algorithm by mclust and regression
#  fit8 <- lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, init_par = "mclust" , seed = 1008)
#  
#  # initialize EM algorithm by random guess
#  fit9 <- lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, init_par = "random" , seed = 1008)

## ---- eval=FALSE--------------------------------------------------------------
#  # summarize LUCID model
#  summary_lucid(fit1)

## ---- eval=FALSE--------------------------------------------------------------
#  # visualze lucid model via a Snakey diagram
#  plot_lucid(fit1)

## ----out.width="70%", echo=FALSE----------------------------------------------
knitr::include_graphics("sankey.png")

## ---- eval=FALSE--------------------------------------------------------------
#  # change node color
#  plot_lucid(fit1, G_color = "yellow")
#  plot_lucid(fit1, Z_color = "red")
#  
#  # change link color
#  plot_lucid(fit1, pos_link_color = "red", neg_link_color = "green")

## ---- eval=FALSE--------------------------------------------------------------
#  # use LUCID model to conduct integrated variable selection
#  # select exposure
#  fit10 <- lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal", K = 2, seed = 1008, Rho_G = 0.1)
#  fit11 <- lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal", K = 2, seed = 1008, Rho_G = seq(0.01, 0.1, by = 0.01))
#  
#  # select omics data
#  fit12 <- lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal", K = 2, seed = 1008, Rho_Z_Mu = 90, Rho_Z_Cov = 0.1, init_par = "random")
#  fit13 <- lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal", K = 2, seed = 1008, Rho_Z_Mu = seq(10, 50, by = 10), Rho_Z_Cov = 0.5, init_par = "random", verbose_tune = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  # tune lucid over a grid of K (note this function may take time to run)
#  tune_lucid <- lucid(G = G, Z = Z, Y = Y_normal, K =2:5)

## ---- eval=FALSE--------------------------------------------------------------
#  # conduct bootstrap resampling
#  boot1 <- boot_lucid(G = G, Z = Z, Y = Y_normal, model = fit1, R = 100)
#  
#  # use 90% CI
#  boot2 <- boot_lucid(G = G, Z = Z, Y = Y_normal, model = fit1, R = 100, conf = 0.9)

## ---- eval=FALSE--------------------------------------------------------------
#  # check distribution for bootstrap replicates of the variable of interest
#  plot(boot1$bootstrap, 1)

## ---- eval=FALSE--------------------------------------------------------------
#  # fit LUCID model with block-wise missing pattern in omics data
#  Z_miss_1 <- Z
#  Z_miss_1[sample(1:nrow(Z), 0.3 * nrow(Z)), ] <- NA
#  fit14 <- lucid(G = G, Z = Z_miss_1, Y = Y_normal, family = "normal", K = 2)
#  
#  # fit LUCID model with sporadic missing pattern in omics data
#  Z_miss_2 <- Z
#  index <- arrayInd(sample(length(Z_miss_2), 0.3 * length(Z_miss_2)), dim(Z_miss_2))
#  Z_miss_2[index] <- NA
#  fit15 <- lucid(G = G, Z = Z_miss_2, Y = Y_normal, family = "normal", K = 2, seed = 1008)
#  
#  # check the imputed omics dataset
#  head(fit15$Z)

