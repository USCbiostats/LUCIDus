# use HELIX Data Challenge data as an example in LUCIDus package
# 2022-02-24
rm(list = ls())



library(tidyverse)
devtools::load_all()

load("~/Documents/Research/LUCID_package/helix_data/exposome.rdata")
load("~/Documents/Research/LUCID_package/helix_data/proteome.rdata")


# 1. Prepare data ====
# exposure and covariates
expo_names = as.character(codebook$variable_name[codebook$family == "Organochlorines"])
cova_names = c("h_mbmi_None", "e3_sex_None", "h_age_None", "h_cohort", "h_edumc_None")
expo = as_tibble(exposome[, c("ID", expo_names)])
cova = as_tibble(covariates[, c("ID", cova_names)])
phen = as_tibble(phenotype[, c("ID", "hs_bmi_c_cat", "hs_zbmi_who")])

# proteomics data
expr_proteome_raw = exprs(proteome)
expr_proteome_sample_id = colnames(expr_proteome_raw)
expr_proteome_feature_id = rownames(expr_proteome_raw)
expr_proteome_annot = fData(proteome)
expr_proteome = as_tibble(t(expr_proteome_raw))
colnames(expr_proteome) = expr_proteome_feature_id
expr_proteome$ID = as.integer(expr_proteome_sample_id)
expr_proteome = expr_proteome[, c(ncol(expr_proteome), 1:(ncol(expr_proteome) - 1))]

# merge data
full_dat = inner_join(expo, cova, by = "ID") %>%
  inner_join(., phen, by = "ID") %>% 
  inner_join(., expr_proteome, by = "ID")
dim(full_dat)

exposure <- full_dat[2:19]
covariate <- full_dat[20:24]
outcome <- full_dat[25:26]
outcome$hs_bmi_c_cat <- ifelse(as.numeric(outcome$hs_bmi_c_cat) < 3,
                               0, 1)
table(outcome$hs_bmi_c_cat)
proteomics <- full_dat[27:62]


fit1 <- est.lucid(G = exposure,
                  Z = scale(proteomics),
                  Y = outcome$hs_zbmi_who,
                  CoY = covariate$h_mbmi_None,
                  family = "normal",
                  K = 2,
                  init_par = "random")
summary(fit1)
pred1 <- predict(fit1,
                 G = exposure,
                 Z = scale(proteomics),
                 Y = outcome$hs_zbmi_who,
                 CoY = covariate$h_mbmi_None)
pred1$pred.y

fit2 <- est.lucid(G = exposure,
                  Z = scale(proteomics),
                  Y = outcome$hs_bmi_c_cat,
                  CoY = covariate$h_mbmi_None,
                  family = "binary",
                  K = 2,
                  init_par = "mclust")
summary(fit2)
pred2 <- predict(fit2,
                 G = exposure,
                 Z = scale(proteomics))




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

tune_mod <- readRDS("../helix_data/tune_mod.rds")
helix_data <- readRDS("~/Documents/Research/LUCID_package/helix_data/helix_data.rds")
usethis::use_data(helix_data)
