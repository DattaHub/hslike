
setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
source("eval_hslike.R")
source("eval_HS.R")

## Example 1
n = 200
A = 7
q = 0.05
theta = c(rep(A,n*q),rep(0,n*(1-q)))
Y = rnorm(n,theta,1)
hsl.fit = eval_hslike(Y)
hs.fit = eval_HS(Y)

library(ggplot2)
est.data = rbind(data.frame(mean = hsl.fit$BetaHat, method = "HS-like",type="Estimate"),
                 data.frame(mean = hs.fit$BetaHat, method = "HS",type="Estimate"),
                 data.frame(mean = theta, method = "True mean",type="Model"),
                 data.frame(mean = Y, method = "Observation",type="Model"))
est.data = cbind(est.data,index=c(1:n))

est.plot = ggplot(est.data,aes(x=index, y=mean, group=method, colour = method)) +
  #geom_line(aes(colour=method),size=1) +
  geom_line(aes(colour=method),position="identity",size=1,show.legend = TRUE)+
  ylab(expression(theta))+
  facet_grid(~type,scale="free")

est.plot <- est.plot+theme_bw()+theme(strip.text.y = element_text(size=12, face="bold")) #+
  #theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

est.plot <- est.plot + theme(axis.title.y = element_text(size = rel(1.5), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.5)))
est.plot <- est.plot + theme(legend.position = "right",legend.key = element_rect(colour = "black"),legend.text=element_text(size=rel(1.2)))
est.plot<- est.plot+ theme(axis.text = element_text(size = rel(1.5)))
print(est.plot)

cairo_ps(file='compare_theta_0119.eps',width=7, height=5)
est.plot
dev.off()

mse1 = sum((hsl.fit$BetaHat-Y)^2)/n
mse2 = sum((hs.fit$BetaHat-Y)^2)/n

mse1 = sum((hsl.fit$BetaHat-hs.fit$BetaHat)^2)/n

## Separate plots
# plot.data = rbind(data.frame(class="Group 1",type="Parameters",value = , x=obswindow),
#                   data.frame(class="Group 2",type="Parameters",value =lambda4(obswindow), x=obswindow),
#                   data.frame(class="Group 1",type="Observations",value =y[,1], x=obswindow),
#                   data.frame(class="Group 2",type="Observations",value =y[,2], x=obswindow),
#                   data.frame(class="Shrinkage",type="Posterior Probabilty",value = post.prob, x = obswindow)
#                   #data.frame(class="Multiscale",type="Posterior Probabilty",value = full_postprob, x = midpts)
# )
# 
# library(ggplot2)
# var.plot <- ggplot(data=plot.data,aes(x=x,y=value,group=class,colour=class))+
#   geom_line(aes(colour=class),position="identity",alpha=0.8,size=1,show.legend = TRUE)+
#   #stat_function(fun=lambda3,colour="red")+stat_function(fun=lambda4,colour="blue")+
#   facet_wrap(~type,ncol=1,scale="free")+theme_bw()+
#   xlab("index")#+
# 
# #scale_colour_manual(values = c("#1b9e77","#d95f02","#7570b3","#e7298a"))
# var.plot <- var.plot + theme(axis.title.y = element_text(size = rel(1), angle = 90))+
#   theme(axis.title.x = element_blank()) #theme(axis.title.x = element_text(size = rel(1)))
# var.plot<- var.plot+ theme(axis.text = element_text(size = rel(1))) +
#   theme(legend.title=element_text(size=12, face="bold"),legend.text=element_text(size=12),
#         legend.position="bottom")
# #theme()
# #theme(legend.title=element_text(size=12, face="bold"),legend.text=element_text(size=12))
# var.plot <- var.plot+theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=10, face="bold"))
# print(var.plot)