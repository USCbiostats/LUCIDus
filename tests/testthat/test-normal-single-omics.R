
test_that("check estimations of LUCID with normal outcome (K = 2)", {
  # run LUCID model
  G <- sim_data$G[1:200, ]
  Z <- sim_data$Z[1:200, ]
  Y_normal <- sim_data$Y_normal[1:200, ]
  cov <- sim_data$Covariate[1:200, ]
  X <- sim_data$X[1:200]
  i <- sample(1:2000, 1)
  cat(paste("seed =", i))
  fit1 <- est.lucid(G = G,
                    Z = Z,
                    Y = Y_normal,
                    CoY = cov,
                    family = "normal",
                    K = 2,
                    seed = i,
                    modelName = "VVV")
  pars <- fit1$pars
  beta_causal <- mean(pars$beta[2, 2:5])
  beta_non <- mean(pars$beta[2, 6:10])
  mu_causal <- mean(abs(pars$mu[1, 1:5] - pars$mu[2, 1:5]))
  mu_non <- mean(abs(pars$mu[1, 6:10] - pars$mu[2, 6:10]))
  gamma_causal <- as.numeric(abs(pars$gamma$beta[1] - pars$gamma$beta[2]))
  gamma_non <- as.numeric(mean(pars$gamma$beta[3:4]))
  sigma <- mean(pars$gamma$sigma)
  
  # check parameters
  expect_equal(beta_causal, log(2), tolerance = 0.2)
  expect_equal(beta_non, 0, tolerance = 0.1)
  expect_equal(mu_causal, 2, tolerance = 0.1)
  expect_equal(mu_non, 0, tolerance = 0.1)
  expect_equal(gamma_causal, 1, tolerance = 0.05)
  expect_equal(gamma_non, 0, tolerance = 0.05)
  expect_equal(sigma, 1, tolerance = 0.05)
})