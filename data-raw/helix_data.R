# use HELIX Data from ISGlobal data challenge as an example in LUCIDus package
# update - 03-18-2022
rm(list = ls())


library(tidyverse)
library(ggplot2)
library(naniar)
devtools::load_all()


# 1. Prepare data ====
load("~/Documents/Research/LUCID_single_omics/exposome_data_challenge/data/exposome.rdata") # exposure without missingness
load("~/Documents/Research/LUCID_single_omics/exposome_data_challenge/data/metabol_serum.rdata")
load("~/Documents/Research/LUCID_single_omics/exposome_data_challenge/data/metabol_urine.rdata")
load("~/Documents/Research/LUCID_single_omics/exposome_data_challenge/data/proteome.rdata")

# 1. individual level data
expo_names = as.character(codebook$variable_name[codebook$family == "Organochlorines"])
cova_names = c("h_mbmi_None", "e3_sex_None", "h_age_None", "h_cohort", "h_edumc_None")
expo = as_tibble(exposome[, c("ID", expo_names)])
cova = as_tibble(covariates[, c("ID", cova_names)])
phen = as_tibble(phenotype[, c("ID", "hs_bmi_c_cat", "hs_zbmi_who")])
# 4. proteome data (1170 x 36 features)
expr_proteome_raw = exprs(proteome)
expr_proteome_sample_id = colnames(expr_proteome_raw)
expr_proteome_feature_id = rownames(expr_proteome_raw)
expr_proteome_annot = fData(proteome)
expr_proteome = as_tibble(t(expr_proteome_raw))
colnames(expr_proteome) = expr_proteome_feature_id
expr_proteome$ID = as.integer(expr_proteome_sample_id)
expr_proteome = expr_proteome[, c(ncol(expr_proteome), 1:(ncol(expr_proteome) - 1))]

full_dat = inner_join(expo, cova, by = "ID") %>%
  inner_join(., phen, by = "ID") %>% 
  inner_join(., expr_proteome, by = "ID")

# select 100 obs
set.seed(123)
index <- sample(1:nrow(full_dat), 100)
dat_select <- full_dat[index, ]
exposure <- dat_select[, c(2:7, 18:19)]
outcome <- dat_select[, 25:26]
outcome$hs_bmi_c_cat <- ifelse(as.numeric(outcome$hs_bmi_c_cat) <= 2,
                               0, 1)
omics <- as_tibble(scale(dat_select[, c(32, 33, 37, 56:62)]))
covariate <- dat_select[, c(20:22)]

helix_data <- list(exposure = exposure,
                   outcome = outcome,
                   omics = omics,
                   covariate = covariate)

usethis::use_data(helix_data, overwrite = TRUE)




# 2. Testing ====
rm(list = ls())
library(LUCIDus)
data("helix_data")
exposome <- helix_data$exposure
proteomics <- helix_data$omics
outcome_norm <- helix_data$outcome["hs_zbmi_who"]
outcome_binary <- helix_data$outcome["hs_bmi_c_cat"]
covariate <- model.matrix(~., helix_data$covariate)[, -1]

## 1 - fit lucid model ====
fit1 <- est.lucid(G = exposome,
                  Z = proteomics,
                  Y = outcome_norm,
                  K = 2)
summary_lucid(fit1)
table(fit1$post.p[, 1] > 0.5)
fit2 <- est.lucid(G = exposome,
                  Z = proteomics,
                  Y = outcome_binary,
                  family = "binary",
                  K = 2)
summary_lucid(fit2)

# unsupervised lucid
fit3 <- est.lucid(G = exposome,
                  Z = proteomics,
                  Y = outcome_norm,
                  K = 2,
                  useY = FALSE)
summary_lucid(fit2)

# include covariates

fit4 <- est.lucid(G = exposome,
                  Z = proteomics,
                  Y = outcome_norm,
                  K = 2,
                  CoY = covariate)
fit5 <- est.lucid(G = exposome,
                  Z = proteomics,
                  Y = outcome_norm,
                  K = 2,
                  CoG = covariate)



# visualize lucid model
plot_lucid(fit1)
plot_lucid(fit1, G_color = "blue")



## 2 - variable selection
fit6 <- est.lucid(G = exposome,
                  Z = proteomics,
                  Y = outcome_norm,
                  K = 2,
                  Rho_G = 0.05)
fit6$select$selectG

fit7 <- est.lucid(G = exposome,
                  Z = proteomics,
                  Y = outcome_norm,
                  K = 2,
                  Rho_Z_Mu = 7,
                  Rho_Z_Cov = 0.1)
fit7$select$selectZ

# refit lucid model with selected variables
fit8 <- est.lucid(G = exposome[, fit6$select$selectG],
                  Z = proteomics[, fit7$select$selectZ],
                  Y = outcome_norm,
                  K = 2)



## 2 - model selection ====
tune_K <- lucid(G = exposome,
                Z = proteomics,
                Y = outcome_norm,
                K = 2:9,
                modelName = NULL)

ggplot(data = tune_K$tune_list,
       aes(x = K, y = BIC, label = round(BIC, 0))) + 
  geom_point() +
  geom_line() + 
  geom_text(hjust = -0.2)
  
# optimal K = 2

# tune optimal penalty in respetive to G
tune_penalty_G <- lucid(G = exposome,
                        Z = proteomics,
                        Y = outcome_norm,
                        K = 2,
                        Rho_G = seq(0.01, 0.05, by = 0.005))
tune_penalty_G$tune_list
tune_result_G <- tune_penalty_G$tune_list
best_penalty_G <- tune_result_G[which.min(tune_result_G$BIC), "Rho_G"]

summary_lucid(tune_penalty_G$best_model)
exposure_select <- tune_penalty_G$best_model$select$selectG

tune_penalty_Z <- lucid(G = exposome,
                        Z = proteomics,
                        Y = outcome_norm,
                        K = 2,
                        Rho_Z_Mu = seq(1, 10, by = 1),
                        Rho_Z_Cov = seq(0.1, 0.5, by = 0.1))
tune_result <- tune_penalty_Z$tune_list
best_Z_penalty <- tune_result[which.min(tune_result$BIC), c("Rho_Z_Mu", "Rho_Z_Cov")]
best_Z_penalty
summary_lucid(tune_penalty_Z$best_model)
omics_select <- tune_penalty_Z$best_model$select$selectZ

# refit the model with selected exposures and omics 
fit5 <- est.lucid(G = exposure[, exposure_select],
                  Z = omics[, omics_select],
                  Y = outcome_norm,
                  K = 2)
summary_lucid(fit5)


## 3 - lucid with missing data ====
# 1 - listwise missing pattern
omics_miss_1 <- omics
omics_miss_1[sample(1:nrow(omics), 0.3 * nrow(omics)), ] <- NA
fit6 <- est.lucid(G = exposure,
                  Z = omics_miss_1,
                  Y = outcome_norm,
                  K = 2)

# 2 - sporadic missing pattern
omics_miss_2 <- as.matrix(omics)
index <- arrayInd(sample(1000, 0.1 * 1000), dim(omics_miss_2))
omics_miss_2[index] <- NA
fit7 <- est.lucid(G = exposure,
                  Z = omics_miss_2,
                  Y = outcome_norm,
                  K = 2)

# 3 - a combination of 2 types of missing pattern
omics_miss_3 <- omics_miss_2
omics_miss_3[sample(1:nrow(omics), 0.2 * nrow(omics)), ] <- NA
vis_miss(as.data.frame(omics_miss_3))

fit8 <- est.lucid(G = exposure,
                  Z = omics_miss_3,
                  Y = outcome_norm,
                  K = 2)

set.seed(1)
proteomics_miss <- as.matrix(proteomics)
index <- arrayInd(sample(1000, 0.1 * 1000), dim(proteomics))
proteomics_miss[index] <- NA # sporadic missing pattern
proteomics_miss[sample(1:100, 20), ] <- NA # listwise missing pattern
vis_miss(as.data.frame(proteomics_miss))

fit6 <- est.lucid(G = exposome,
                  Z = proteomics_miss,
                  Y = outcome_norm,
                  K = 2)
fit7 <- est.lucid(G = exposome,
                  Z = proteomics_miss,
                  Y = outcome_binary,
                  K = 2)





## 4 - inference ====
set.seed(123)
boot1 <- boot.lucid(G = exposome, 
                    Z = proteomics,
                    Y = outcome_norm,
                    model = fit1,
                    R = 100)
plot(boot1$bootstrap, 10)


summary_lucid(fit1, boot.se = boot1)


## 5 - prediction ====
pred1 <- predict_lucid(G = exposome,
                       Z = proteomics,
                       Y = outcome_norm,
                       model = fit1)
pred2 <- predict_lucid(G = exposome,
                       Z = proteomics,
                       # Y = outcome_norm,
                       model = fit1)
table(pred1$pred.x)
table(pred2$pred.x)




