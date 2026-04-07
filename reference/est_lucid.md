# Fit LUCID model to conduct integrated clustering

The Latent Unknown Clustering with Integrated Data (LUCID) performs
integrative clustering using multi-view data. LUCID model is estimated
via EM algorithm for model-based clustering. It also features variable
selection, integrated imputation, bootstrap inference and visualization
via Sankey diagram.

## Usage

``` r
est_lucid(
  G,
  Z,
  Y,
  CoG = NULL,
  CoY = NULL,
  K = 2,
  family = c("normal", "binary"),
  useY = TRUE,
  tol = 0.001,
  max_itr = 1000,
  max_tot.itr = 10000,
  Rho_G = 0,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  modelName = NULL,
  seed = 123,
  init_impute = c("mclust", "lod"),
  init_par = c("mclust", "random"),
  verbose = FALSE
)
```

## Arguments

- G:

  Exposures, a numeric vector, matrix, or data frame. Categorical
  variable should be transformed into dummy variables. If a matrix or
  data frame, rows represent observations and columns correspond to
  variables.

- Z:

  Omics data, a numeric matrix or data frame. Rows correspond to
  observations and columns correspond to variables.

- Y:

  Outcome, a numeric vector. Categorical variable is not allowed. Binary
  outcome should be coded as 0 and 1.

- CoG:

  Optional, covariates to be adjusted for estimating the latent cluster.
  A numeric vector, matrix or data frame. Categorical variable should be
  transformed into dummy variables.

- CoY:

  Optional, covariates to be adjusted for estimating the association
  between latent cluster and the outcome. A numeric vector, matrix or
  data frame. Categorical variable should be transformed into dummy
  variables.

- K:

  Number of latent clusters. An integer greater or equal to 2. User can
  use [`lucid`](lucid.md) to determine the optimal number of latent
  clusters.

- family:

  Distribution of outcome. For continuous outcome, use "normal"; for
  binary outcome, use "binary". Default is "normal".

- useY:

  Flag to include information of outcome when estimating the latent
  cluster. Default is TRUE.

- tol:

  Tolerance for convergence of EM algorithm. Default is 1e-3.

- max_itr:

  Max number of iterations for EM algorithm.

- max_tot.itr:

  Max number of total iterations for `est_lucid` function. `est_lucid`
  may conduct EM algorithm for multiple times if the algorithm fails to
  converge.

- Rho_G:

  A scalar. This parameter is the LASSO penalty to regularize exposures.
  If user wants to tune the penalty, use the wrapper function `lucid`

- Rho_Z_Mu:

  A scalar. This parameter is the LASSO penalty to regularize
  cluster-specific means for omics data (Z). If user wants to tune the
  penalty, use the wrapper function `lucid`

- Rho_Z_Cov:

  A scalar. This parameter is the graphical LASSO penalty to estimate
  sparse cluster-specific variance-covariance matrices for omics data
  (Z). If user wants to tune the penalty, use the wrapper function
  `lucid`

- modelName:

  The variance-covariance structure for omics data. See
  [`mclust::mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.html)
  for details.

- seed:

  An integer to initialize the EM algorithm or imputing missing values.
  Default is 123.

- init_impute:

  Method to initialize the imputation of missing values in LUCID.
  "mclust" will use `mclust:imputeData` to implement EM Algorithm for
  Unrestricted General Location Model to impute the missing values in
  omics data; `lod` will initialize the imputation via relacing missing
  values by LOD / sqrt(2). LOD is determined by the minimum of each
  variable in omics data.

- init_par:

  Method to initialize the EM algorithm. "mclust" will use mclust model
  to initialize parameters; "random" initialize parameters from uniform
  distribution.

- verbose:

  A flag indicates whether detailed information for each iteration of EM
  algorithm is printed in console. Default is FALSE.

## Value

A list which contains the several features of LUCID, including:

- pars:

  Estimates of parameters of LUCID, including beta (effect of exposure),
  mu (cluster-specific mean for omics data), sigma (cluster-specific
  variance-covariance matrix for omics data) and gamma (effect estimate
  of association between latent cluster and outcome)

- K:

  Number of latent cluster

- modelName:

  Geometric model to estiamte variance-covariance matrix for omics data

- likelihood:

  The log likelihood of the LUCID model

- post.p:

  Posterior inclusion probability (PIP) for assigning observation i to
  latent cluster j

- Z:

  If missing values are observed, this is the complet dataset for omics
  data with missing values imputed by LUCID

## References

Cheng Peng, Jun Wang, Isaac Asante, Stan Louie, Ran Jin, Lida Chatzi,
Graham Casey, Duncan C Thomas, David V Conti, A Latent Unknown
Clustering Integrating Multi-Omics Data (LUCID) with Phenotypic Traits,
Bioinformatics, btz667, https://doi.org/10.1093/bioinformatics/btz667.

## Examples

``` r
if (FALSE) { # \dontrun{
# use simulated data
G <- sim_data$G
Z <- sim_data$Z
Y_normal <- sim_data$Y_normal
Y_binary <- sim_data$Y_binary
cov <- sim_data$Covariate

# fit LUCID model with continuous outcome
fit1 <- est_lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, 
seed = 1008)

# fit LUCID model with block-wise missing pattern in omics data
Z_miss_1 <- Z
Z_miss_1[sample(1:nrow(Z), 0.3 * nrow(Z)), ] <- NA
fit2 <- est_lucid(G = G, Z = Z_miss_1, Y = Y_normal, family = "normal", K = 2)

# fit LUCID model with sporadic missing pattern in omics data
Z_miss_2 <- Z
index <- arrayInd(sample(length(Z_miss_2), 0.3 * length(Z_miss_2)), dim(Z_miss_2))
Z_miss_2[index] <- NA
# initialize imputation by imputing 
fit3 <- est_lucid(G = G, Z = Z_miss_2, Y = Y_normal, family = "normal", 
K = 2, seed = 1008, init_impute = "lod") 
LOD
# initialize imputation by mclust
fit4 <- est_lucid(G = G, Z = Z_miss_2, Y = Y, family = "normal", K = 2, 
seed = 123, init_impute = "mclust") 

# fit LUCID model with binary outcome
fit5 <- est_lucid(G = G, Z = Z, Y = Y_binary, family = "binary", K = 2,
seed = 1008)

# fit LUCID model with covariates
fit6 <- est_lucid(G = G, Z = Z, Y = Y_binary, CoY = cov, family = "binary", 
K = 2, seed = 1008)

# use LUCID model to conduct integrated variable selection
# select exposure
fit6 <- est_lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal", 
K = 2, seed = 1008, Rho_G = 0.1)
# select omics data
fit7 <- est_lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal",
K = 2, seed = 1008, Rho_Z_Mu = 90, Rho_Z_Cov = 0.1, init_par = "random")

} # }
```
