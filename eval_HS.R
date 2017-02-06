eval_HS <- function(Y)
{
  TauTrue <- 0.1
  p <- length(Y)
  n <- 1
  Beta = rep(0,p)
  Sigma2 = 1
  Sigma = 1
  Tau = TauTrue
  Lambda = rep(1,p)
  Sigma = sqrt(Sigma2)
  nmc = 10000
  burn = 100
  BetaSave = matrix(0, nrow=nmc, ncol=p)
  LambdaSave = matrix(0, nrow=nmc, ncol=p)
  TauSave = rep(0, nmc)
  Sigma2Save = rep(0, nmc)
  Res2 = Y
  for(t in 1:(nmc+burn))
  {
    if(t %% 1000 == 0) cat("Iteration ",t, "\n")
    # First block-update Beta
    a = (Tau^2)*(Lambda^2)
    b = n*a
    s = sqrt(Sigma2*a/{1+b})
    m = {b/{1+b}}*Y
    Beta = rnorm(p, m, s)
    Theta = Beta/(Sigma*Lambda)
    # Now update Sigma2
    # Jeffreys prior is assumed
    Res2 = {Y^2}/{1+(Tau^2)*{Lambda^2}}
    RSS = sum(Res2)
    Sigma2 = 1/rgamma(1,n*p/2, rate = RSS/2)
    Sigma = sqrt(Sigma2)
    # Now update Tau^2 using slice sampling
    eta = 1/(Tau^2)
    u = runif(1,0,1/(eta + 1))
    ub = (1-u)/u
    a = (p+1)/2
    b = sum(Theta^2)/2
    ub2 = pgamma(ub,a,rate=b)
    u2 = runif(1,0,ub2)
    eta = qgamma(u2,a,rate=b)
    Tau = 1/sqrt(eta)
    # Now update Lambda, comment out for global l shrinkage only
    Z = Y/(Sigma*Theta)
    V2 = 1/rgamma(p,1,rate=(Lambda^2+1)/2)
    num1 = n*V2*(Theta^2)
    den = 1 + num1
    s = sqrt(V2/den)
    m = {num1/den}*Z
    Lambda = rnorm(p,m,s)
    if(t > burn)
    {
      BetaSave[t-burn,] = Beta
      LambdaSave[t-burn,] = Lambda
      TauSave[t-burn] = Tau
      Sigma2Save[t-burn] = Sigma2
    }
  }
  BetaHat = apply(BetaSave,2,mean)
  LambdaHat = apply(abs(LambdaSave),2,mean)
  TauHat = mean(TauSave)
  Sigma2Hat = mean(Sigma2Save)
  return(list(BetaHat=BetaHat,Lambda=LambdaHat,Tau=TauHat,Sigma=Sigma2Hat))
}

  