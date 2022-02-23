# use HELIX Data Challenge data as an example in LUCIDus package
# 2022-02-22
rm(list = ls())
devtools::load_all()

dat <- readRDS("../helix_data/helix_data.rds")
# prepare data
exposure <- dat$exposure
omics <- scale(dat$omics)
outcome_norm <- dat$outcome[, 1]
covariate <- model.matrix(~., dat$covariate)[, -1]
tune_mod <- lucid(G = exposure,
                  Z = omics,
                  Y = outcome_norm,
                  CoY = covariate,
                  K = 2:4,
                  Rho_Z_Mu = c(30, 40, 50),
                  Rho_Z_Cov = c(0.1, 0.2, 0.3),
                  Rho_G = c(0.005, 0.01),
                  seed = 123)
tune_mod$tune_list
summary(tune_mod$best_model)

saveRDS(tune_mod,
        file = "../helix_data/tune_mod.rds")
