# Deprecated function est.lucid

This function deprecates. Please use est_lucid instead.

## Usage

``` r
est.lucid(
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
  modelName = "VVV",
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
