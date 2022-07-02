# test plot_lucid function

test_that("check whether plot_lucid function could work", {
  # run LUCID model
  G <- sim_data$G[1:200, ]
  Z <- sim_data$Z[1:200, ]
  Y_normal <- sim_data$Y_normal[1:200, ]
  cov <- sim_data$Covariate[1:200, ]
  X <- sim_data$X[1:200]
  i <- 1008
  invisible(capture.output(fit1 <- est_lucid(G = G,
                                             Z = Z,
                                             Y = Y_normal,
                                             CoY = cov,
                                             family = "normal",
                                             K = 2,
                                             seed = i,
                                             modelName = "VVV")))
  plot1 <- plot_lucid(fit1)
  plot2 <- plot_lucid(fit1, G_color = "black", X_color = "red", 
                      pos_link_color = "green", fontsize = 10)
  expect_equal(class(plot1), c("sankeyNetwork", "htmlwidget"))
  expect_equal(class(plot2), c("sankeyNetwork", "htmlwidget"))
})