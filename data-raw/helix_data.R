# use HELIX Data Challenge data as an example in LUCIDus package
# 2022-02-21
library(Biobase)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mclust)
library(reshape2)
library(LUCIDus)
library(R.utils) # keep track the time
library(networkD3)

# 1. Prepare data ====
load(file = "../helix_data/exposome.rdata")
load(file = "../helix_data/metabol_serum.rdata")
load(file = "../helix_data/metabol_urine.rdata")
load(file = "../helix_data/proteome.rdata")


# 1-1. extract exposure, outcome and covaraites
expo_names <- as.character(codebook$variable_name[codebook$family == "Organochlorines"])
cova_names <- c("h_mbmi_None", "e3_sex_None", "h_age_None", "h_cohort", "h_edumc_None")
expo <- exposome[, c("ID", expo_names)]
cova <- covariates[, c("ID", cova_names)]
phen <- phenotype[, c("ID", "hs_bmi_c_cat", "hs_zbmi_who")]

# 1-2 extract omics data
# serum data (1198 x 177 features)
expr_serum_raw = exprs(metabol_serum)
expr_serum_sample_id = colnames(expr_serum_raw)
expr_serum_feature_id = paste0("serum_", rownames(expr_serum_raw))
expr_serum = as_tibble(t(expr_serum_raw))
colnames(expr_serum) = expr_serum_feature_id
expr_serum$ID = as.integer(expr_serum_sample_id)
expr_serum = expr_serum[, c(ncol(expr_serum), 1:(ncol(expr_serum) - 1))]
# urine data (1192 x 44 features): 
expr_urine_raw = exprs(metabol_urine)
expr_urine_sample_id = colnames(expr_urine_raw)
expr_urine_feature_id = paste0("urine_", rownames(expr_urine_raw))
expr_urine = as_tibble(t(expr_urine_raw)) 
colnames(expr_urine) = expr_urine_feature_id
expr_urine$ID = as.integer(expr_urine_sample_id)
expr_urine = expr_urine[, c(ncol(expr_urine), 1:(ncol(expr_urine) - 1))]
# proteome data (1170 x 36 features)
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
  inner_join(., expr_serum, by = "ID") %>% 
  inner_join(., expr_urine, by = "ID") %>% 
  inner_join(., expr_proteome, by = "ID")
dim(full_dat) # 1152 obs x 283 variables (ID, 36 + 44 + 177 = 257 metabolites, 2 pheno, 5 cova, 18 exposures)
sum(is.na(full_dat)) # 0 NA value

# write.csv(full_dat, 
#           file = "../helix_data/dat_organochlorines.csv", 
#           quote = FALSE, 
#           row.names = FALSE)
exposure <- full_dat[, 2:19]
outcome <- full_dat[, 25:26]
covariate <- full_dat[, 20:24]
serum_metab <- full_dat[, 27:203]
urine_metab <- full_dat[, 204:247]
protein <- full_dat[, 248:283]



