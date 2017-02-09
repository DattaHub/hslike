source("https://bioconductor.org/biocLite.R")
biocLite("multtest")
require(multtest)
data(golub)
class(golub)
dim(golub)
head(golub)

## Gene names
dim(golub.gnames)
golub.gnames[1:4, ]
## Class labels
golub.cl

teststat = mt.teststat(golub, golub.cl,test="t.equalvar")
teststat = qnorm(pt(teststat,df=36))
#teststat = (teststat - mean(teststat))/sd(teststat)

qqnorm(teststat)
qqline(teststat)

require(ggplot2)
plt = ggplot(data.frame(teststat), aes(sample = teststat)) + stat_qq() + theme_bw()
plt

########### 

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

y = teststat
n = length(teststat)

## Fit MCMC
ptm<-proc.time()
evalHS=eval_HS(y)
tMCMC<-proc.time() - ptm
coefHS=evalHS$BetaHat
tau_HS=evalHS$Tau

kappa.hs <- 1/(1+ evalHS$Lambda^2*evalHS$Tau^2)
coefHS <- ifelse(kappa.hs<0.5,evalHS$BetaHat,0)

## Fit using Horseshoe package - buggy code, doesn't converge !!
# library(horseshoe)
# tau.example <- HS.MMLE(y, var(y))
# coefHS<-HS.post.mean(y, tau.example, 1)

## Fit using my implementation !
# source("logshrink.R")
# hsl_fit <- eval_hslike(data,nmc=10000)
# coefHS<-colMeans(hsl_fit$ThetaSave)

##Fit EM
ptm<-proc.time()
evalem=eval_EM(y)
tEM<-proc.time() - ptm
coefEM = evalem$thetahat
a_EM = evalem$a
#nzeroem=length(which(coefEM==0 & theta==0))


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

#------------------#
## Plots by JD ##

est.data <- rbind(data.frame(mean=coefEM,method="HS-Mode",x=y),
                  #data.frame(mean=coefHS,method="HS-Mean",x=y),
                  data.frame(mean=coefSCAD,method="SCAD",x=y),
                  data.frame(mean=coefMCP,method="MCP",x=y),
                  data.frame(mean=coefLASSO,method="LASSO",x=y))

est.plot = ggplot(est.data,aes(x=x, y=mean, group=method, colour = method)) +
  #geom_line(aes(colour=method),size=1) +
  geom_line(aes(colour=method),position="identity",size=1,show.legend = TRUE)+
  ylab(expression(theta))+xlab("Test Statistics")
  #facet_grid(~type,scale="free")
est.plot <- est.plot+theme_bw()+theme(strip.text.y = element_text(size=12, face="bold")) #+
#theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

est.plot <- est.plot + theme(axis.title.y = element_text(size = rel(1.5), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.5)))
est.plot <- est.plot + theme(legend.position = c(0.9,0.4),legend.key = element_rect(colour = "black"),legend.text=element_text(size=rel(1.2)))
est.plot<- est.plot+ theme(axis.text = element_text(size = rel(1.5)))
print(est.plot)

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
cairo_ps(file='compare_thetahat_leuk_2.eps',width=7, height=5)
est.plot
dev.off()

## Multiple testing 

rawp = 2 * (1 - pnorm(abs(teststat)))
adjusted = mt.rawp2adjp(rawp, c("BH","Bonferroni","BY"))
apply(mt.reject(adjusted$adj[order(adjusted$index), ], 0.05)$which,2,sum)
##--- Non-zeroes---##
# myf <- function(x) { length(x[x>0]) }
# # mutate(est.data, value = myf(mean))
# apply(est.data,1,function(x){ifelse(is.numeric(x),length(x[x>0]),NA)})
# ##
length(coefEM[coefEM==0])
length(coefHS[coefHS==0])
length(coefSCAD[coefSCAD==0])
length(coefMCP[coefMCP==0])
length(coefLASSO[coefLASSO==0])

adjusted$adjp[1:10, ]
adjusted$adj[order(adjusted$index)[1:10], ]

which = mt.reject(adjusted$adj[order(adjusted$index), ], 0.05)$which[, 2]
#golub.gnames[which, 2]
head(golub.gnames[which, 2])

### Posterior mean and multiple testing plot 

est.data2 <- cbind(est.data, result = factor(which,labels=c("Not rejected","Rejected")))
est.plot<- ggplot(est.data2,aes(x=x, y=mean, group=method)) + 
  #geom_line(aes(colour=method),size=1) +
  geom_line(aes(colour=result),position="identity",size=1,show.legend = TRUE)+
  ylab(expression(theta))+xlab("Test Statistics")+
  facet_wrap(~method,nrow=2)

est.plot <- est.plot+theme_bw()+theme(strip.text.y = element_text(size=12, face="bold")) +
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

est.plot <- est.plot + theme(axis.title.y = element_text(size = rel(1.5), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.5)))
est.plot <- est.plot + theme(legend.position = "bottom",legend.key = element_rect(colour = "black"),legend.text=element_text(size=rel(1.2)),legend.title=element_text(size=rel(1.5)))
est.plot<- est.plot+ theme(axis.text = element_text(size = rel(1.5)))
print(est.plot)

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
cairo_ps(file='facet_compare_thetahat_leuk_2.eps',width=7, height=5)
est.plot
dev.off()

## Supplementary plot

teststat = mt.teststat(golub, golub.cl,test="t.equalvar")
teststat = qnorm(pt(teststat,df=36))
#teststat = (teststat - mean(teststat))/sd(teststat)
df = data.frame(TS = teststat)

n = length(teststat)
plt = ggplot(df, aes(x = df$TS)) + 
  geom_histogram(aes(y = ..density..),binwidth=0.15,colour="black",fill="gray") +
  stat_function(fun=dnorm,size=rel(1.2),linetype="dashed")+
  stat_function(fun=dnorm, args = list(mean=mean(df$TS),sd = sd(df$TS)))+
  xlab("Test statistics") # +
  # labs(title="Histogram of z-values and N(0,1) curve") 

plt <- plt+theme_bw()+theme(strip.text.y = element_text(size=12, face="bold")) #+
#theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

plt <- plt + theme(axis.title.y = element_text(size = rel(1.5), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.5)))
plt <- plt + theme(legend.position = "right",legend.key = element_rect(colour = "black"),legend.text=element_text(size=rel(1.2)))
plt<- plt+ theme(axis.text = element_text(size = rel(1.5)))
print(plt)

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
cairo_ps(file='compare_density_leuk_0205.eps',width=7, height=5)
plt
dev.off()

### Another way 

n        <- length(teststat)
mean     <- mean(teststat)
sd       <- sd(teststat)
binwidth <- 0.25

qplot(teststat, geom = "histogram", breaks = seq(-5, 5, binwidth), 
      colour = I("black"), fill = I("white"),
      xlab = "Test statistics", ylab = "Count") +
  # Create normal curve, adjusting for number of observations and binwidth
  stat_function( 
    fun = function(x, mean, sd, n, bw){ 
      dnorm(x = x, mean = mean, sd = sd) * n * bw
    }, 
    args = c(mean = 0, sd = 1, n = n, bw = binwidth))


