# use HELIX Data from ISGlobal data challenge as an example in LUCIDus package
# update - 03-18-2022
rm(list = ls())


library(tidyverse)
library(Biobase)
devtools::load_all()


helix_data <- readRDS("../helix_data/helix_data.rds")
helix_data$exposure <- helix_data$exposure[, c(1:6, 17, 18)]
helix_data$omics <- scale(helix_data$omics[, 36:45])
helix_data$outcome$hs_bmi_c_cat <- ifelse(as.numeric(helix_data$outcome$hs_bmi_c_cat) <= 2,
                                          0, 1)
table(helix_data$outcome$hs_bmi_c_cat)
usethis::use_data(helix_data, overwrite = TRUE)


covariate <- model.matrix(~., helix_data$covariate)[, -1]

# 2. Testing ====
## 1. fit lucid model ====
fit1 <- est.lucid(G = exposure,
                  Z = omics,
                  Y = outcome_norm,
                  K = 2)
summary_lucid(fit1)
fit2 <- est.lucid(G = exposure,
                  Z = omics,
                  Y = outcome_norm,
                  K = 2,
                  useY = FALSE)
summary_lucid(fit2)

# include covariates
fit3 <- est.lucid(G = exposure,
                  Z = omics,
                  Y = outcome_norm,
                  K = 2,
                  CoY = covariate)
summary_lucid(fit3)

fit4 <- est.lucid(G = exposure,
                  Z = omics,
                  Y = outcome_norm,
                  K = 2,
                  CoG = covariate)
summary_lucid(fit4)


# visualize lucid model
plot_lucid(fit1)
plot_lucid(fit1, G_color = "blue")


## 2 - model selection ====
tune_K <- lucid(G = exposure,
                Z = omics,
                Y = outcome_norm,
                K = 2:6)
plot(x = tune_K$tune_list$K, 
     y = tune_K$tune_list$BIC, 
     type = "b", 
     xlab = "K", 
     ylab = "BIC")
# optimal 

tune_penalty_G <- lucid(G = helix_data$exposure,
                        Z = omics,
                        Y = outcome_norm,
                        K = 2,
                        Rho_G = c(0.001, 0.003, 0.005, 0.007, 0.01))

tune_penalty_G$tune_list
summary_lucid(tune_penalty_G$best_model)

# refit the model with selected exposures
exposure_select <- tune_penalty_G$best_model$select$selectG
fit_5 <- est.lucid(G = helix_data$exposure[, exposure_select],
                   Z = omics,
                   Y = outcome_norm,
                   K = 2)
summary_lucid(fit_5)
table(fit_final$post.p[, 1] > 0.5)

saveRDS(tune_mod,
        file = "../helix_data/tune_mod.rds")




