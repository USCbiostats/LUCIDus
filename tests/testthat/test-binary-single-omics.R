# LUCID - single omics, binary outcome


test_that("check estimations of LUCID with binary outcome (K = 2)", {
  # run LUCID model
  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y_binary <- sim_data$Y_binary[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  i <- 1008
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- est.lucid(G = G,
                                             Z = Z,
                                             Y = Y_binary,
                                             CoY = cov,
                                             family = "binary",
                                             K = 2,
                                             seed = i,
                                             useY = TRUE,
                                             modelName = "VVV")))
  pars <- fit1$pars
  beta_causal <- mean(pars$beta[2, 2:5])
  beta_non <- mean(pars$beta[2, 6:10])
  mu_causal <- mean(abs(pars$mu[1, 1:5] - pars$mu[2, 1:5]))
  mu_non <- mean(abs(pars$mu[1, 6:10] - pars$mu[2, 6:10]))
  gamma <- as.numeric(pars$gamma$beta)
  
  # check parameters
  expect_equal(beta_causal, log(2), tolerance = 0.2)
  expect_equal(beta_non, 0, tolerance = 0.1)
  expect_equal(mu_causal, 2, tolerance = 0.1)
  expect_equal(mu_non, 0, tolerance = 0.1)
  expect_equal(gamma, c(-0.5, 0.9, 0.8, -0.8), tolerance = 0.2)
})


test_that("check variable selection on G", {
  # run LUCID model
  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y_binary <- sim_data$Y_binary[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  # cat(paste("test2 - seed =", i, "\n"))
  i <- 1008
  invisible(capture.output(fit1 <- est.lucid(G = G,
                                             Z = Z,
                                             Y = Y_binary,
                                             CoY = cov,
                                             family = "binary",
                                             K = 2,
                                             seed = i,
                                             useY = TRUE,
                                             modelName = "VVV",
                                             Rho_G = 0.1)))
  
  # check parameters
  expect_equal(class(fit1$select$selectG), "logical")
  expect_equal(as.vector(fit1$select$selectG), 
               c(rep(TRUE, 4), rep(FALSE, 6)))
})


test_that("check variable selection on Z", {
  # run LUCID model
  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y_binary <- sim_data$Y_binary[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  # cat(paste("test3 - seed =", i, "\n"))
  i <- 1008
  invisible(capture.output(fit1 <- est.lucid(G = G,
                                             Z = Z,
                                             Y = Y_binary,
                                             CoY = cov,
                                             family = "binary",
                                             K = 2,
                                             seed = i,
                                             useY = TRUE,
                                             modelName = "VVV",
                                             Rho_Z_Mu =  25,
                                             Rho_Z_Cov = 0.15)))
  
  # check parameters
  expect_equal(class(fit1$select$selectG), "logical")
})


test_that("check whether arguments of est.lucid work", {
  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y_binary <- sim_data$Y_binary[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  # cat(paste("test4 - seed =", i, "\n"))
  i <- 1008
  invisible(capture.output(fit1 <- est.lucid(G = G,
                                             Z = Z,
                                             Y = Y_binary,
                                             CoY = cov,
                                             family = "binary",
                                             K = 2,
                                             seed = i,
                                             useY = TRUE,
                                             modelName = NULL)))
  invisible(capture.output(fit2 <- est.lucid(G = G,
                                             Z = Z,
                                             Y = Y_binary,
                                             CoY = cov,
                                             family = "binary",
                                             K = 2,
                                             seed = i,
                                             useY = TRUE,
                                             init_par = "random")))
  expect_equal(class(fit1$modelName), "character")
  expect_equal(fit2$init_par, "random")
})


test_that("check whether est.lucid throws an erorr with continuous outcome or outcome coded rather than 0 and 1", {
  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y_binary <- sim_data$Y_binary[1:500, ]
  Y_normal <- sim_data$Y_normal[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  i <- 1008
  # cat(paste("test4 - seed =", i, "\n"))
  expect_error(est.lucid(G = G,
                         Z = Z,
                         Y = Y_normal,
                         CoY = cov,
                         family = "binary",
                         K = 2,
                         seed = i,
                         useY = TRUE,
                         modelName = NULL))
  expect_error(est.lucid(G = G,
                         Z = Z,
                         Y = Y_binary + 1,
                         CoY = cov,
                         family = "binary",
                         K = 2,
                         seed = i,
                         useY = TRUE,
                         modelName = NULL))
})