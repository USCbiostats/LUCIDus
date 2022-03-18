# use HELIX Data from ISGlobal data challenge as an example in LUCIDus package
# update - 03-18-2022
rm(list = ls())


library(tidyverse)
library(Biobase)
devtools::load_all()

load("~/Documents/Research/LUCID_package/helix_data/exposome.rdata")
load("~/Documents/Research/LUCID_package/helix_data/proteome.rdata")


# 1. Prepare data ====
# # exposure and covariates
# expo_names = as.character(codebook$variable_name[codebook$family == "Organochlorines"])
# cova_names = c("h_mbmi_None", "e3_sex_None", "h_age_None", "h_cohort", "h_edumc_None")
# expo = as_tibble(exposome[, c("ID", expo_names)])
# cova = as_tibble(covariates[, c("ID", cova_names)])
# phen = as_tibble(phenotype[, c("ID", "hs_bmi_c_cat", "hs_zbmi_who")])
# 
# # proteomics data
# expr_proteome_raw = exprs(proteome)
# expr_proteome_sample_id = colnames(expr_proteome_raw)
# expr_proteome_feature_id = rownames(expr_proteome_raw)
# expr_proteome_annot = fData(proteome)
# expr_proteome = as_tibble(t(expr_proteome_raw))
# colnames(expr_proteome) = expr_proteome_feature_id
# expr_proteome$ID = as.integer(expr_proteome_sample_id)
# expr_proteome = expr_proteome[, c(ncol(expr_proteome), 1:(ncol(expr_proteome) - 1))]
# 
# # merge data
# full_dat = inner_join(expo, cova, by = "ID") %>%
#   inner_join(., phen, by = "ID") %>%
#   inner_join(., expr_proteome, by = "ID")
# dim(full_dat)
# 
# 
# # prepare data
# exposure <- full_dat[2:19]
# covariate <- full_dat[20:24]
# outcome <- full_dat[25:26]
# outcome$hs_bmi_c_cat <- ifelse(as.numeric(outcome$hs_bmi_c_cat) < 3,
#                                0, 1)
# table(outcome$hs_bmi_c_cat)
# proteomics <- full_dat[27:62]
# helix_data <- list(exposure = exposure,
#                    outcome = outcome,
#                    proteomics = scale(proteomics),
#                    covariate = covariate)


helix_data <- readRDS("../helix_data/helix_data.rds")
exposure <- helix_data$exposure
omics <- scale(helix_data$omics[, c(36:45)])
outcome_norm <- helix_data$outcome[, 1]
covariate <- model.matrix(~., helix_data$covariate)[, -1]

# 2. Testing ====
## 1. fit lucid model ====
fit1 <- est.lucid(G = exposure,
                  Z = omics,
                  Y = outcome_norm,
                  K = 4)
summary_lucid(fit1)
fit2 <- est.lucid(G = exposure,
                  Z = omics,
                  Y = outcome_norm,
                  K = 4,
                  useY = FALSE,
                  verbose = TRUE)
summary_lucid(fit2)




## 2 - model selection ====
tune_mod <- lucid(G = exposure,
                  Z = omics,
                  Y = outcome_norm,
                  K = 2:6,
                  # Rho_Z_Mu = c(30, 40, 50),
                  # Rho_Z_Cov = c(0.1, 0.2, 0.3),
                  # Rho_G = c(0.005, 0.01),
                  seed = 123,
                  modelName = "VVV",
                  useY = TRUE)
plot(x = tune_mod$tune_list$K, y = tune_mod$tune_list$BIC, 
     type = "b", xlab = "K", ylab = "BIC")


tune_mod$tune_list
summary_lucid(tune_mod$best_model)
best_tune <- tune_mod$tune_list[which.min(tune_mod$tune_list$BIC), ]
# refit the model with selected variables
selectG <- tune_mod$best_model$select$selectG
selectZ <- tune_mod$best_model$select$selectZ
exposure_select <- exposure[, selectG]
omics_select <- omics[, selectZ]

fit_final <- est.lucid(G = exposure_select,
                       Z = omics_select,
                       Y = outcome_norm,
                       CoY = covariate,
                       K = 2,
                       useY = TRUE)
summary_lucid(fit_final)
table(fit_final$post.p[, 1] > 0.5)

saveRDS(tune_mod,
        file = "../helix_data/tune_mod.rds")

# tune_mod <- readRDS("../helix_data/tune_mod.rds")
# helix_data <- readRDS("../helix_data/helix_data.rds")
usethis::use_data(helix_data, overwrite = TRUE)

