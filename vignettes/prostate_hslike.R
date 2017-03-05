#setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
load(prostz.RDa)

load("prostz.Rda")
y = prostz
cat(length(y),"\n") # correct data should have 6,033 entries

teststat = y
qqnorm(teststat)
qqline(teststat)

require(ggplot2)
plt = ggplot(data.frame(teststat), aes(sample = teststat)) + stat_qq() + theme_bw()
plt

##
source("eval_HS.R")
source("eval_EM.R")

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

est.data <- rbind(#data.frame(mean=coefEM,method="HS-Mode",x=y),
                  data.frame(mean=coefHS,method="HS-Mean",x=y),
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
est.plot <- est.plot + theme(legend.position = c(0.2,0.8),legend.key = element_rect(colour = "black"),legend.text=element_text(size=rel(1.2)))
est.plot<- est.plot+ theme(axis.text = element_text(size = rel(1.5)))
print(est.plot)

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
cairo_ps(file='compare_thetahat_prost.eps',width=7, height=5)
est.plot
dev.off()

## Multiple testing 
library(multtest)
rawp = 2 * (1 - pnorm(abs(teststat)))
adjusted = mt.rawp2adjp(rawp, c("Bonferroni", "Holm", "Hochberg", "SidakSS", 
                                "SidakSD", "BH", "BY", "ABH", "TSBH"))
apply(mt.reject(adjusted$adj[order(adjusted$index), ], 0.05)$which,2,sum)
##--- Non-zeroes---##
# myf <- function(x) { length(x[x>0]) }
# # mutate(est.data, value = myf(mean))
# apply(est.data,1,function(x){ifelse(is.numeric(x),length(x[x>0]),NA)})
# ##
n-length(coefEM[coefEM==0])
n-length(coefHS[coefHS==0])
n-length(coefSCAD[coefSCAD==0])
n-length(coefMCP[coefMCP==0])
n-length(coefLASSO[coefLASSO==0])

adjusted$adjp[1:10, ]
adjusted$adj[order(adjusted$index)[1:10], ]

which = mt.reject(adjusted$adj[order(adjusted$index), ], 0.05)$which[, 7]
#golub.gnames[which, 2]
#head(golub.gnames[which, 2])

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

#setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
#cairo_ps(file='facet_compare_thetahat_prost.eps',width=7, height=5)
# est.plot
# dev.off()

## Supplementary plot

# teststat = mt.teststat(golub, golub.cl,test="t.equalvar")
# teststat = qnorm(pt(teststat,df=36))
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
cairo_ps(file='compare_density_prost.eps',width=7, height=5)
plt
dev.off()