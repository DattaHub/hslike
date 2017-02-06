library(pracma)
library(Bessel)
#setwd("C:/Users/jyotishka/Documents/R/HSPlus")
setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HSPlus")
high <- 3
low <- -3
nsample <- 1e5
nplot <- 1e3

# Horseshoe
# lamb <- abs(rcauchy(nsample))
# hs_smpl <- rnorm(nsample,mean = 0, sd = sqrt(lamb))
# hs_den <- density(hs_smpl,from=low,to=high,n=nsample)

theta_seq = seq(low,high,len = nplot)
hsval <- rep(0,nplot)
robval <- rep(0,nplot)
for (i in 1:nplot){
  t = theta_seq[i]
  prior <- function(x) dnorm(t,mean=0,sd=x)*2*dcauchy(x,location=0, scale=1) 
  hsv <- integrate(prior,lower=0,upper=Inf)
  hsval[i] <- hsv$value
  #   robust <- function(x) dnorm(t,mean=0,sd=sqrt(x))*densigamma(x,alpha=0.5,beta=0.5) 
  robust <- function(x) dnorm(t,mean=0,sd=x)*4/(pi)^2*log(x)/(x^2-1)
  rbv <- integrate(robust,0,1000)
  robval[i]<- rbv$value
}

#cauchy

tval <- dt(theta_seq,df = 1)
#gdp
a = 1
b = 1
gdpval <- (1/2*b)*(1+abs(theta_seq)/(a*b))^(-(a+1))
#laplace
laplval <- (1/2)*exp(-abs(theta_seq))
#DL 
a = 0.5
DLval <- 1/(2^(0.5+a/2)*gamma(a))*(abs(theta_seq))^((a-1)/2)*besselK(sqrt(2*abs(theta_seq)),1-a)

tau = 1; a = 2*tau^2
HLval <- 1/(2*pi*sqrt(a))*log(1+a/theta_seq^2)

eps = 3/4
HL <- function(x) ifelse(x>=1,abs(x)^(1-eps)*log(1+2/x^2),abs(x)^(-(1-eps))*log(1+2/x^2))
C1 <- integrate(HL,Inf,0)
C2 <- integrate(HL,0,Inf)
HLval = HL(theta_seq)/(C1$value+C2$value)

ntype <- 7
colors <- c("blue","red","darkgreen","black","magenta","cyan","brown")
linetype <- c(1:ntype) 
plotchar <- seq(16,16+ntype,1)

plot(gdpval~theta_seq,type = 'l',lwd=2,
     lty=linetype[1], col=colors[1], pch=plotchar[1],xlim=c(low,high),
     ylim=c(0,1),ylab=expression(pi(theta)),xlab=expression(theta)) 
lines(theta_seq,tval, lwd=2,
      lty=linetype[2], col=colors[2], pch=plotchar[2]) 
lines(theta_seq,laplval, lwd=2,
      lty=linetype[3], col=colors[3], pch=plotchar[3]) 
lines(theta_seq,hsval, lwd=2,
      lty=linetype[4], col=colors[4], pch=plotchar[4]) 
lines(theta_seq,robval, lwd=2,
      lty=linetype[5], col=colors[5], pch=plotchar[5])
lines(theta_seq,DLval, lwd=2,
      lty=linetype[6], col=colors[6], pch=plotchar[6])
lines(theta_seq,HLval, lwd=2,
      lty=linetype[7], col=colors[7], pch=plotchar[7])
legend(2,0.9, c("GDP","Cauchy","Laplace","HS","HSPlus","DL","HSLike"), cex=0.8, col=colors,
       lty=linetype,xjust=0.5)

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
dev.copy2pdf(file = "dens_zero_all_0129.pdf")
dev.copy2eps(file = "dens_zero_all_0129.eps")