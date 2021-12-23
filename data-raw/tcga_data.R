# codes to prepare the TCGA data ====

# use the lusc data from package whitening
# install.packages("whitening)
library(whitening)
data("lusc")

## 1. create DNA methylation data ====
# choose the top 10 DNA methylation sites that are most significantly associated
# with survtime
dim(lusc$methyl)
methyl_name <- colnames(lusc$methyl)
dat <- as.data.frame(cbind(lusc$survivalTime, lusc$methyl))
colnames(dat)[1] <- "survtime"
p_res <- rep(0, dim(lusc$methyl)[2])
for(i in 1:dim(lusc$methyl)[2]) {
  temp_fit <- lm(as.formula(paste0("survtime~", methyl_name[i])), data = dat)
  p_res[i] <- summary(temp_fit)$coefficients[2, 4]
}
summary(p_res)
methyl_name_top10 <- methyl_name[rank(p_res) <= 10]
methyl_name_bottom10 <- methyl_name[rank(p_res) > (234 - 10)]
Z <- as.data.frame(lusc$methyl)[c(methyl_name_top10, methyl_name_bottom10)]


## 2. impute missing data in exposure (smoke pack) ====
G <- lusc$packs
sum(is.na(G)) # 20 observations are missing
# assign random values to smoking packs
set.seed(1008)
G[is.na(G)] <- sample(G[!is.na(G)], 20)


# test the lusc data
# G <- as.matrix(G)
# Z <- as.matrix(Z)
# Y <- as.matrix(lusc$survivalTime)
# fit1 <- est.lucid(G = G,
#                   Z = Z,
#                   Y = Y,
#                   family = "normal",
#                   K = 2,
#                   seed = 1027)
# summary(fit1)
survtime <- as.matrix(lusc$survivalTime)
colnames(survtime) <- "survtime"
lusc_data <- list(smoke_pack = G,
                  methyl = Z,
                  survtime = survtime)
usethis::use_data(lusc_data, overwrite = TRUE)
