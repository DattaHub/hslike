##Horseshoe and approximate horseshoe
## A. Bhadra, J. Datta, N. G. Polson, B. Willard
## November 14, 2016
## Function description:
## eval_HS(y): Full horseshoe MCMC. Sample (tau, lambda_i | y_i). Sigma = 1. Global-local shrinkage.
## eval_OHagan(y): Horseshoe MCMC with fixed tau=1, Sigma=1. Sample (lambda_i | y_i). Pure local shrinkage.
## eval_EM(y): Evaluate posterior mode by "approximate" penalty (upper and lower bounded by a=4, a=2). tau=fixed at 1

rm(list=ls())
set.seed(100)

##Full horseshoe MCMC (sample tau, lambda_i given y)
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
  nmc = 15000
  burn = 5000
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
    # Res2 = {Y^2}/{1+(Tau^2)*{Lambda^2}}
    # RSS = sum(Res2)
    # Sigma2 = 1/rgamma(1,n*p/2, rate = RSS/2)
    # Sigma = sqrt(Sigma2)
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

####eval_E-M: Evaluate posterior mode. 
eval_EM<-function(Y)
{
	a<-5 #a=2 upper bound, a=4 lower bound
	p<-length(Y)
	tol = 1e-10
	err=1
	t=1
	thetahat = y
	while (err>tol){
		#thetahatnew = y/(1 + (1/(pi*(a^1.5)))*(a/(thetahat^2) - (a/(thetahat^2 +a)) ))
		#anew = (1/(p*pi*sqrt(a)))*sum(a-a*thetahatnew^2/(a+thetahatnew^2))
		thetahatnew = y/(1 + (sqrt(a))/(pi*thetahat^2*(thetahat^2+a)))
		anew = (a^1.5/(p*pi))*sum(1/(thetahatnew^2+a))
		err<-sum((thetahatnew-thetahat)^2) + (a-anew)^2
		cat("EM Iteration ",t, "Error ", err, "a ", a, "\n")
		thetahat=thetahatnew
		a=anew
		t<-t+1
	}
	return(list(thetahat=thetahat, a=a))
}


##sparse theta
n=1000
theta= c(rep(3,10),rep(-3,10), rep(0,n-20))
#theta=rep(10,n)
y=rnorm(n=n,mean=theta,sd=1)

## Fit MCMC
ptm<-proc.time()
evalHS=eval_HS(y)
tMCMC<-proc.time() - ptm
coefHS=evalHS$BetaHat
tau_HS=evalHS$Tau

##Fit EM
ptm<-proc.time()
evalem=eval_EM(y)
tEM<-proc.time() - ptm
coefEM = evalem$thetahat
a_EM = evalem$a
nzeroem=length(which(coefEM==0 & theta==0))


##Fit LASSO
library(glmnet)
ptm<-proc.time()
LASSOfit=cv.glmnet(x=diag(n), y=y)
tLASSO<-proc.time() - ptm
coefLASSO=coef(LASSOfit,s="lambda.min")
coefLASSO=coefLASSO[-1]

##Fit MCP
library(sparsenet)
ptm<-proc.time()
MCPfit=cv.sparsenet(x=diag(n), y=y)
tMCP<-proc.time() - ptm
coefMCP=coef(MCPfit)
coefMCP=coefMCP[-1]

##Fit SCAD
library(ncvreg)
ptm<-proc.time()
SCADfit=cv.ncvreg(X=diag(n), y=y, penalty="SCAD")
tSCAD<-proc.time() - ptm
coefSCAD=coef(SCADfit)
coefSCAD=coefSCAD[-1]


# loglik_HS = sum(dnorm(y,mean=coefHS,sd=1,log=T))
# #loglik_OHagan = sum(dnorm(y,mean=thetafit_OHagan,sd=1,log=T))
# loglik_EM = sum(dnorm(y,mean=coefEM,sd=1,log=T))
# #loglik_EMa = sum(dnorm(y,mean=coefEMa,sd=1,log=T))
# loglik_LASSO = sum(dnorm(y,mean=coefLASSO,sd=1,log=T))
# loglik_MCP = sum(dnorm(y,mean=coefMCP,sd=1,log=T))
# loglik_SCAD = sum(dnorm(y,mean=coefSCAD,sd=1,log=T))
# loglik_true=sum(dnorm(y,mean=theta,sd=1,log=T))

mse_HS=sum((theta-coefHS)^2)
mse_EM=sum((theta-coefEM)^2)
mse_LASSO=sum((theta-coefLASSO)^2)
mse_MCP=sum((theta-coefMCP)^2)
mse_SCAD=sum((theta-coefSCAD)^2)

cat("\n RESULTS: \n")
cat("===========\n")
# #Print MSE
cat("\n \t \t \t & mode  & mean  & scad  & mcp  & lasso \\\\ \n", "MSE \t \t \t&", round(mse_EM,2), "&", round(mse_HS,2), "&", round(mse_SCAD,2), "&", round(mse_MCP,2), "&", round(mse_LASSO,2), "\\\\ \n ")

## #Print subset selection results: detection of zeros
EM0=length(which(coefEM==0 & theta==0))
MCP0=length(which(coefMCP==0 &theta==0))
SCAD0=length(as.vector(which(coefSCAD==0 &theta==0)))
LASSO0=length(which(coefLASSO==0 & theta==0))
cat("Cor\\_Z \t \t",  "&", EM0, "& NA &", SCAD0, "&", MCP0, "&", LASSO0, "\\\\  \n")

## #Print subset selection results: detection of non-zeros
EM1=length(which(coefEM!=0 & theta!=0))
MCP1=length(which(coefMCP!=0 &theta!=0))
SCAD1=length(as.vector(which(coefSCAD!=0 &theta!=0)))
LASSO1=length(which(coefLASSO!=0 & theta!=0))
cat("Cor\\_Nz \t", " &", EM1, "& NA &", SCAD1, "&", MCP1, "&", LASSO1, "\\\\  \n")

# # Print times
cat("Time \t \t  \t &", tEM[[1]], "&", tMCMC[[1]], "&", tSCAD[[1]], "&", tMCP[[1]], "&", tLASSO[[1]], "\\\\  \n")



## # Plots
par(mfrow=c(1,5))
par(mar=c(5.1,5.5,4.1,2.1))
par(pch=16)
plot(theta, coefEM, col="blue", xlab=expression(theta), ylab=expression(hat(theta)), xlim=c(-6,6), ylim = c(-6,6), cex.lab=2, cex.axis =1.5)
points(theta, coefHS, col="red")
points(theta, coefSCAD, col="green")
points(theta, coefMCP, col="cyan")
points(theta, coefLASSO, col="black")
abline(b=1, a=0, lwd=1.5, lty=2)
title(expression(paste("Plot of true ", theta, " vs. ", hat(theta))), cex.main=3)
legend("bottomright", c("HS_MODE", "HS_MEAN", "SCAD", "MCP", "LASSO"), pch=c(16,16,16,16,16), col=c("blue", "red", "green","cyan", "black"), cex=1.5)


plot(y,coefEM, col="blue",  xlab="Y", ylab=expression(hat(theta)), xlim=c(-6,6), ylim = c(-6,6), cex.lab=2, cex.axis =1.5)
points(y, coefHS, col="red")
points(y, coefSCAD, col="green")
points(y,coefMCP, col="cyan")
points(y,coefLASSO, col="black")
title(expression(paste("Plot Y vs. ", hat(theta))), cex.main=3)
abline(b=1, a=0, lwd=1.5, lty=2)
legend("bottomright", c("HS_MODE", "HS_MEAN", "SCAD", "MCP", "LASSO"), pch=c(16,16,16,16,16), col=c("blue", "red", "green","cyan", "black"), cex=1.5)


plot(coefEM, coefSCAD, xlab=expression(hat(theta)[MODE]), ylab= expression(hat(theta)[SCAD]), cex.lab=2, cex.axis=1.5 )
title(expression(paste("Plot of ", hat(theta)[MODE], " vs. ", hat(theta)[SCAD])), cex.main=3)

plot(coefEM, coefMCP, xlab=expression(hat(theta)[MODE]), ylab= expression(hat(theta)[MCP]), cex.lab=2, cex.axis=1.5 )
title(expression(paste("Plot of ", hat(theta)[MODE], " vs. ", hat(theta)[MCP])), cex.main=3)

plot(coefEM, coefLASSO, xlab=expression(hat(theta)[MODE]), ylab= expression(hat(theta)[LASSO]), cex.lab=2, cex.axis=1.5)
title(expression(paste("Plot of ", hat(theta)[MODE], " vs. ", hat(theta)[LASSO])), cex.main=3)



  