##Horseshoe and approximate horseshoe
## A. Bhadra, J. Datta, N. G. Polson, B. Willard
## January 25, 2017
## Function description: eval_EM_regression(y,x) - Evaluate posterior mode by "approximate" penalty.

rm(list=ls())
set.seed(100)

####eval_E-M: Evaluate posterior mode for theta in the model: Y |X, theta = N(X*theta, 1); p(theta)=horseshoe;
### Arguments: Y = n X 1 vector of observations; X = n X p matrix of predictors
### Result: thetahat= p X 1 vector of regression coefficients; ahat = posterior mode for global parameter
eval_EM_regression<-function(Y,X)
{
	a<-0.05
	p<-ncol(X)
	n<-length(Y)
	tol = 1e-8
	err=1
	t=1
	thetahat = rep(1,p)
	thetahat=scale(t(X)%*%Y)
	while (err>tol){
		utilde = (1/(2*pi))*(a^1.5/(thetahat^2*(thetahat^2 +a)))
		D=diag(as.vector(a/(2*utilde)))
		DXT = D%*%t(X)
		XTY = t(X)%*%Y
		cov=D - DXT%*%(solve(X%*%DXT + diag(n)))%*%t(DXT)
		thetahatnew = cov%*%XTY
		anew = (a^1.5/(p*pi))*sum(1/(thetahatnew^2+a))
		err<-sum((thetahatnew-thetahat)^2) + (a-anew)^2
        cat("EM Iteration ",t, "Error ", err, "a ", a, "\n")
		thetahat=thetahatnew
		a=anew
		t<-t+1
	}
	return(list(thetahat=thetahat, ahat=a))
}

## Compare posterior mode with lasso and sparsenet
## Set parse theta and generate data
n= 70
p= 350 #n=60, p=300 good
theta= c(rep(3,10),rep(-3,10), rep(0,p-20))

## Case 1: X matrix is i.i.d standard normal (uncorrelated X)
x=matrix(rnorm(n*p), nrow=n, ncol=p)
# ## #Case 2: X matrix comes from a "k" factor model (correlated X)
# k <- 1
# B <- matrix(rep(0.1,p*k),nrow=p,ncol=k)
# x <- matrix(0,n,p)
  # for (i in (1:n)) {
    # Factor <- matrix(rnorm(k,mean=0,sd=1))
    # x[i,] <- B%*%Factor+rnorm(p,mean=0,sd=0.1)
# }

## Generate y
y=rnorm(n=n,mean=x%*%theta,sd=1)

## Run EM
ptm<-proc.time()
evalem=eval_EM_regression(Y=y,X=x)
coefEM = evalem$thetahat
a_EM = evalem$ahat
tEM<-proc.time() - ptm
nzeroem=length(which(coefEM==0 & theta==0))
Yhat_EM=x%*%coefEM

##Fit LASSO
library(glmnet)
ptm<-proc.time()
LASSOfit=cv.glmnet(x=x, y=y)
tLASSO<-proc.time() - ptm
coefLASSO=coef(LASSOfit,s="lambda.min")
coefLASSO=coefLASSO[-1]
Yhat_LASSO=x%*%coefLASSO

##Fit MCP
library(sparsenet)
ptm<-proc.time()
MCPfit=cv.sparsenet(x=x, y=y)
tMCP<-proc.time() - ptm
coefMCP=coef(MCPfit)
coefMCP=coefMCP[-1]
Yhat_MCP=x%*%coefMCP

##Fit SCAD
library(ncvreg)
ptm<-proc.time()
SCADfit=cv.ncvreg(X=x, y=y, penalty="SCAD")
tSCAD<-proc.time() - ptm
coefSCAD=coef(SCADfit)
coefSCAD=coefSCAD[-1]
Yhat_SCAD=x%*%coefSCAD

## Run HS Postmean
library(horseshoe)
ptm<-proc.time()
evalHS=horseshoe(y=y,X=x, method.tau="halfCauchy",burn=2000,nmc=5000)
coefHS = evalHS$BetaHat
tMCMC<-proc.time() - ptm
nzeroHS=length(which(coefHS==0 & theta==0))
Yhat_HS=x%*%coefHS

##MSE of the estimates
mse_EM=sum((theta-coefEM)^2)
mse_LASSO=sum((theta-coefLASSO)^2)
mse_MCP=sum((theta-coefMCP)^2)
mse_SCAD=sum((theta-coefSCAD)^2)
mse_HS = sum((theta-coefHS)^2)
 
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

## Time
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


plot(y,Yhat_EM, col="blue",  xlab="Y", ylab=expression(hat(Y)), xlim=c(-6,6), ylim = c(-6,6), cex.lab=2, cex.axis =1.5)
points(y, Yhat_HS, col="red")
points(y, Yhat_SCAD, col="green")
points(y,Yhat_MCP, col="cyan")
points(y,Yhat_LASSO, col="black")
title(expression(paste("Plot Y vs. ", hat(Y))), cex.main=3)
abline(b=1, a=0, lwd=1.5, lty=2)
legend("bottomright", c("HS_MODE", "HS_MEAN", "SCAD", "MCP", "LASSO"), pch=c(16,16,16,16,16), col=c("blue", "red", "green","cyan", "black"), cex=1.5)


plot(coefEM, coefSCAD, xlab=expression(hat(theta)[MODE]), ylab= expression(hat(theta)[SCAD]), cex.lab=2, cex.axis=1.5 )
title(expression(paste("Plot of ", hat(theta)[MODE], " vs. ", hat(theta)[SCAD])), cex.main=3)

plot(coefEM, coefMCP, xlab=expression(hat(theta)[MODE]), ylab= expression(hat(theta)[MCP]), cex.lab=2, cex.axis=1.5 )
title(expression(paste("Plot of ", hat(theta)[MODE], " vs. ", hat(theta)[MCP])), cex.main=3)

plot(coefEM, coefLASSO, xlab=expression(hat(theta)[MODE]), ylab= expression(hat(theta)[LASSO]), cex.lab=2, cex.axis=1.5)
title(expression(paste("Plot of ", hat(theta)[MODE], " vs. ", hat(theta)[LASSO])), cex.main=3)




  