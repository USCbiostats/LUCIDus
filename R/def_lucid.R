# control parameters for EM algorithm
def.control <- function(tol = 1e-4,
                        max_itr = 1000,
                        max_tot.itr = 10000){
  structure(list(tol = tol,
                 max_itr = max_itr,
                 max_tot.itr = max_tot.itr))
}


# define tuning parameters
def.tune <- function(Rho_G = 0,
                     Rho_Z_InvCov = 0,
                     Rho_Z_CovMu = 0,
                     Select_G = FALSE,
                     Select_Z = FALSE) {
  if(Select_G){
    Rho_G == Rho_G
  }
  if(Select_Z){
    Rho_Z_InvCov = Rho_Z_InvCov
    Rho_Z_CovMu = Rho_Z_CovMu
  }
  structure(list(Select_G = Select_G,
                 Select_Z = Select_Z,
                 Rho_G = Rho_G,
                 Rho_Z_InvCov = Rho_Z_InvCov,
                 Rho_Z_CovMu = Rho_Z_CovMu))
}



