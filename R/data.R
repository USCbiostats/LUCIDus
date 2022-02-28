#' @title A simulated dataset for LUCID 
#'
#' @description This is an example dataset to illustrate LUCID model. It is simulated
#' by assuming there are 2 latent clusters in the data. We assume the exposures 
#' are associated with latent cluster which ultimately affects the PFAS concentration
#' and liver injury in children. The latent clusters are also characterized by 
#' differential levels of metabolites.
#' 
#' @format A list with 5 matrices corresponding to exposures (G), omics data (Z),
#' a continuous outcome, a binary outcome and 2 covariates (can be used either 
#' as CoX or CoY). Each matrice contains 2000 observations.
#' \describe{
#'     \item{G}{10 exposures}
#'     \item{Z}{10 metabolites}
#'     \item{Y_normal}{Outcome, PFAS concentration in children}
#'     \item{Y_bninary}{Bianry outcome, liver injury status}
#'     \item{Covariates}{2 continous covariates, can be treated as either CoX or 
#'                       CoY}
#'     \item{X}{Latent clusters}
#' }
"sim_data"




#' @title TCGA LUSC data
#' 
#' @description This is a modified DNA methylation data combining selected 
#' clinical exposure (number of smoking packs per year) and outcome (survival 
#' time) for 130 patients with lung squamous cell carcinoma (LUSC). The original 
#' .RData file is extracted from R package 'whitening' (Strimmer et al. 2021), 
#' which is a pre-processed sample from The Cancer Genome Atlas (TCGA) database.
#' Please note: the missing values in smoking packs are imputed by random filling 
#' so the dataset is served as a illustration of LUCID analysis and shouldn't be
#' used for formal analysis.
#' 
#' @format A list with 3 matrices corresponding to exposure (G), 
#' omics measurements (Z) and a continuous outcome (survival time). 
#' \describe{
#'     \item{smoke_pack}{exposure, number of smoking packs per year}
#'     \item{methyl}{20 methylation measurements}
#'     \item{survtime}{survival time}
#' }
#' 
#' @source https://cran.r-project.org/web/packages/whitening/index.html
"lusc_data"