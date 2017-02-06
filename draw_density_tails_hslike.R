library(pracma)
library(Bessel)
#setwd("C:/Users/jyotishka/Documents/R/HSPlus")
#setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HSPlus")
high <- 25
low <- 10
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

# tau = 1; a = 4*tau^2
# HLval <- 1/(2*pi*sqrt(a))*log(1+a/theta_seq^2)

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
     ylim=c(0,0.006),ylab=expression(pi(theta)),xlab=expression(theta)) 
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
legend(20,0.006, c("GDP","Cauchy","Laplace","HS","HSPlus","DL","HSLike"), cex=0.8, col=colors,
       lty=linetype,xjust=0.5)

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
dev.copy2pdf(file = "dens_tail_all_0129.pdf")
dev.copy2eps(file = "dens_tail_all_0129.eps")



### Local shrinkage prior density
high <- 8
low <- 4
#nsample <- 1e5
nplot <- 1e2

theta_seq = seq(low,high,len = nplot)
islash = (sqrt(2/pi))*(1-exp(-1/(2*theta_seq^2)))
slash = (sqrt(1/pi))*(1-exp(-theta_seq^2))/(theta_seq^2)
tval <- dt(theta_seq,df = 1)
#laplace
laplval <- (1/2)*exp(-abs(theta_seq))

ntype <- 3
colors <- c("blue","red","darkgreen","black","magenta","cyan","brown")
linetype <- c(1:ntype) 
plotchar <- seq(16,16+ntype,1)

plot(islash~theta_seq,type = 'l',lwd=2,
     lty=linetype[1], col=colors[1], pch=plotchar[1],xlim=c(low,high),
     ylim=c(0,0.04),ylab=expression(pi(theta)),xlab=expression(theta)) 
lines(theta_seq,tval, lwd=2,
      lty=linetype[2], col=colors[2], pch=plotchar[2]) 
lines(theta_seq,laplval, lwd=2,
      lty=linetype[3], col=colors[3], pch=plotchar[3]) 
legend(7,0.03, c("Inverse Slash","Cauchy","Laplace"), cex=1.2, col=colors,
       lty=linetype,xjust=0.5)

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
dev.copy2pdf(file = "local_tails_slash.pdf")
dev.copy2eps(file = "local_tails_slash.eps")

origin.priors.data = rbind(data.frame(prior="Cauchy", prob=tval),
                           data.frame(prior="Laplace", prob=laplval),
                           data.frame(prior="iSlash", prob=islash))
origin.priors.data = cbind(origin.priors.data, theta=theta_seq)

# origin.plot = ggplot(origin.priors.data, 
#                      aes(x=theta, y=prob, group=prior)) +
#   geom_line(aes(linetype=prior),size=1) + theme_bw()+
#   ylab(expression(pi(lambda))) + xlab(expression(lambda))
# print(origin.plot)

origin.plot = ggplot(origin.priors.data, 
                     aes(x=theta, y=prob, group=prior)) +
  geom_line(aes(linetype=prior),size = 1)+ylab("")+ ylab(expression(pi(lambda))) + xlab(expression(lambda))

origin.plot <- origin.plot+theme_bw()+theme(legend.position="none")+#scale_y_discrete(breaks="")+
  scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash","twodash"))

origin.plot <- origin.plot+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

origin.plot <- origin.plot + theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.2)))
origin.plot <- origin.plot + theme(legend.position = "right",legend.key = element_rect(colour = "black"),
                                   legend.text=element_text(size=rel(1.2)),legend.title=element_text(size=rel(1.2)),
                                   legend.key.size=unit(1.2,"cm"))
origin.plot<- origin.plot+ theme(axis.text = element_text(size = rel(1.2)))
# labs(title="Comparison of Different Priors")+theme(plot.title = element_text(size = rel(3)))
# labs <- list("Laplace(1)","Cauchy","Inverse-Slash")
# origin.plot<- origin.plot+ scale_linetype_manual(values=1:5,labels=labs)
origin.plot <- origin.plot + theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.2)))
origin.plot<- origin.plot+ theme(axis.text = element_text(size = rel(1.2)))
origin.plot <- origin.plot+theme(strip.text.x = element_text(size=15, face="bold"),
                                 strip.text.y = element_text(size=15, face="bold"))

print(origin.plot)
setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
ggsave("local-tails-2.eps", origin.plot, width=7, height=5)




## Prior on Kappa

high <- 1;low <- 0
nsample <- 1e5;nplot <- 1e3

kappa_seq = seq(low,high,len = nplot)
kapval = sqrt(1/(2*pi))*(1-exp(-kappa_seq/(2*(1-kappa_seq))))*kappa_seq^(-1.5)*(1-kappa_seq)^(-0.5)
betaval = dbeta(kappa_seq,shape1 = 0.5, shape2 = 0.5)


ntype <- 2
colors <- c("blue","red","darkgreen","black","magenta","cyan","brown")
linetype <- c(1:ntype) 
plotchar <- seq(16,16+ntype,1)
plot(kapval~kappa_seq,type = 'l',lwd=2,
     lty=linetype[1], col=colors[1], pch=plotchar[1],xlim=c(low,high),
     ylim = c(0,10), ylab=expression(pi(theta)),xlab=expression(theta)) 
lines(kappa_seq,betaval, lwd=2,lty=linetype[2], col=colors[2], pch=plotchar[2]) 
legend(0.5,8, c("HSLike","Beta(1/2,1/2)"), cex=1.2, col=colors,
       lty=linetype,xjust=0.5)

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
dev.copy2pdf(file = "compare_kappa_hslike.pdf")
dev.copy2eps(file = "compare_kappa_hslike.eps")
