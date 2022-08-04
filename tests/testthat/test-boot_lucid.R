# LUCID - bootstrap

test_that("check functionality of boot_lucid", {
  # run LUCID model
  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y_binary <- sim_data$Y_binary[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  i <- 1008
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- est_lucid(G = G,
                                             Z = Z,
                                             Y = Y_binary,
                                             CoY = cov,
                                             family = "binary",
                                             K = 2,
                                             seed = i,
                                             useY = TRUE,
                                             modelName = "VVV")))
  invisible(capture.output(boot1 <- boot_lucid(G = G,
                                               Z = Z,
                                               Y = Y_binary,
                                               CoY = cov,
                                               model = fit1,
                                               R = 50)))
  expect_equal(class(boot1), "list")
  expect_equal(names(boot1), c("beta", "mu", "gamma", "bootstrap"))
  expect_equal(class(boot1$bootstrap), "boot")
})
