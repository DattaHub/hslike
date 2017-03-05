#'
#'  Slice sampler for the Normal observations model with HS-like prior.
#'  TODO: More detail on priors.
#'
#'  @param Y observations
#'  @return A list of samples
#'  @export
#'
eval_hslike <- function(Y, nmc=100, burn=100, 
                        eff_zero=1e-5, verbose=TRUE) {
  n <- length(Y)
  Sigma2 = 1
  Sigma = sqrt(Sigma2)
  Tau = 1
  Theta = rep(1, n)
  Lambda = rep(1, n)
  Nu = rep(1, n)
  
  ThetaSave = matrix(0, nrow=nmc, ncol=n)
  LambdaSave = matrix(0, nrow=nmc, ncol=n)
  TauSave = rep(0, nmc)
  Sigma2Save = rep(0, nmc)
  
  for(t in 1:(nmc+burn)) {
    if(isTRUE(verbose) && t %% 200 == 0) 
      cat("Iteration ",t, "\n")
    
    lam_div_tau = Lambda / Tau
    
    # First block-update Theta
    # b = 2 * lam_div_tau^2
    # s = sqrt(Sigma2/{1+b})
    # m = Y / (1 + b)
    # Beta = rnorm(n, m, s)
    # Theta_new = Beta * Lambda
    lam_tau_weight = 1. / (1 + lam_div_tau^2)
    Theta_new = rnorm(n, lam_tau_weight * Y, 
                      sqrt(lam_tau_weight * Sigma2)) 
    
    stopifnot(!any(is.nan(Theta_new)))
    Theta = Theta_new
    
    # Now update Sigma2
    # Jeffreys prior is assumed
    # Res2 = {Y^2}/{1+(2*Lambda^2)/(Tau^2)}
    # Res2 = Y^2 / (1 + 2 * lam_div_tau)
    #Res2 = (Y-Beta)^2
    # RSS = sum(Res2)
    ## Should comment?
    # RSS = crossprod(Y - Theta)
    # Sigma2 = 1 / rgamma(1,(n + 1)/2, rate = RSS/2)
    # Sigma = sqrt(Sigma2)
    
    ## Now update Tau^2 using slice sampling
    eta = 1 / Tau^2
    u = runif(1, 0, 1/(eta + 1))
    ub = (1 - u) / u
    a = (n + 1)/ 2
    b = 0.5*crossprod(Lambda * Theta)
    ub2 = pgamma(ub, a, rate=b)
    u2 = runif(1, 0, ub2)
    eta = qgamma(u2, a, rate=b)
    Tau_new = 1/sqrt(eta)
    
    # XXX: A hack to avoid 0/0 in Lambda / Tau.
    Tau_new = ifelse(Tau_new < eff_zero, eff_zero, Tau_new)
    
    # DEBUG: Remove. 
    #cat(sprintf("t=%d, Tau=%g \n", t, Tau))
    stopifnot(!any(is.nan(Tau_new)))
    Tau = Tau_new
    
    # Now update Lambda, comment out for global shrinkage only
    beta_div_tau = Theta / Tau
    Lambda2_new = rgamma(n, shape = 3/2, rate=0.5*beta_div_tau^2 + Nu/2)
    
    # XXX: A hack to avoid 0/0 in Lambda / Tau.
    Lambda2_new = ifelse(Lambda2_new < eff_zero, eff_zero, Lambda2_new)
    
    # DEBUG: Remove. 
    #cat(sprintf("t=%d, mean(Lambda2 > 0 (eff))=%g \n", t, 
    #            mean(Lambda2_new > eff_zero)))
    stopifnot(!any(is.nan(Lambda2_new)))
    Lambda2 = Lambda2_new
    
    Lambda = sqrt(Lambda2)
    
    # Now Update Nu
    # Nu = rgamma(n, shape=3/2, rate=Lambda2/2) 
    #Nu = rexp(n, rate=Lambda2)
    
    Nu = (-2/Lambda2)*log(1-runif(n)*(1-exp(-Lambda2/2)))
    
    #if(anyNA(Nu)){ browser() }
    
    if(t > burn)
    {
      ThetaSave[t-burn,] = Theta
      LambdaSave[t-burn,] = Lambda
      TauSave[t-burn] = Tau
      Sigma2Save[t-burn] = Sigma2
    }
  }
  
  # ThetaHat = apply(ThetaSave, 2, mean)
  # LambdaHat = apply(abs(LambdaSave), 2, mean)
  # TauHat = mean(TauSave)
  # Sigma2Hat = mean(Sigma2Save)
  return(list(ThetaSave=ThetaSave, 
              LambdaSave=LambdaSave, 
              TauSave=TauSave, 
              Sigma2Save=Sigma2Save 
              # ThetaHat=ThetaHat,
              # LambdaHat=LambdaHat, 
              # TauHat=TauHat, 
              # Sigma2Hat=Sigma2Hat
  ))
}