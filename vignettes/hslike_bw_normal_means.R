
setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
source("logshrink.R")
set.seed(198)
theta_1 = c(rep(6, 10), rep(0, 90))
Y_1 = rnorm(100, theta_1, 1)
#Y_2 = rnorm(100, theta_1, 3)

theta_2 = c(rep(7, 10), rep(3,10), rep(0, 80))
Y_3 = rnorm(100, theta_2, 1)
#Y_4 = rnorm(100, theta_2, 3)

prob_set = list(list(name='theta_1, sigma=1', obs=Y_1),
                list(name='theta_2, sigma=1', obs=Y_3))  

theta.data = lapply(prob_set,
                    function(set) 
                      list('settings'=set$name,
                           'samples'=eval_hslike(set$obs, 
                                                 nmc=5000)))

theta.data = lapply(theta.data,
                    function(sample.data)
                      cbind(melt(sample.data$samples$ThetaSave), 
                            'settings'=sample.data$settings, 
                            'var'='theta'))

theta.data = Reduce(rbind, theta.data)

colnames(theta.data)[1:2] = c('sample', 'component')

obs.data = lapply(prob_set, 
                  function(set) 
                    data.frame(component=seq_along(set$obs), 
                               value=set$obs,
                               sample=NA, 
                               settings=set$name, 
                               var="obs")) 
obs.data = Reduce(rbind, obs.data)

all.data = rbind(obs.data, theta.data)

plot.data = all.data %>% 
  group_by(component, settings, var) %>%
  summarise(upper = quantile(value, prob=0.95), 
            lower = quantile(value, prob=0.05), 
            middle = median(value))

theta.plot = ggplot(plot.data,
                    aes(x=component, y=middle, group=component,
                        colour=var)) + theme_bw() +   
  geom_pointrange(aes(ymin=lower, ymax=upper), size=0.2, alpha=0.5) + 
  facet_grid(settings ~ ., scales="free_y") + #, labeller = label_both) +
  scale_colour_manual(values=c("#D55E00","#0072B2")) +
  xlab("") + ylab("")

print(theta.plot)
setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
ggsave("theta_mcmc_runs-2.pdf", plot=theta.plot)


## Posterior mean plot 
#Plot the posterior mean for a range of deterministic values

#Example with 20 signals, rest is noise
#Posterior mean for the signals is plotted in blue
par(mfrow=c(1,2))
#truth <- c(rep(0, 80), rep(8, 20))
truth <- c(rep(0, 80), rep(3,10), rep(6, 10))
data <- truth + rnorm(100)
tau.example <- HS.MMLE(data, 1)
plot(data, HS.post.mean(data, tau.example, 1),
     col = c(rep("black", 80), rep("red",10), rep("blue", 20)),main="HS",ylab="theta")
# hsfit <- eval_HS(data)
# post.mean.theta = hsfit$BetaHat/sqrt(hsfit$Sigma*hsfit$Lambda^2)
# plot(data, post.mean.theta,col = c(rep("black", 80), rep("blue", 20)))
hsl_fit <- eval_hslike(data,nmc=5000)
plot(data, colMeans(hsl_fit$ThetaSave),col = c(rep("black", 80), rep("red",10), rep("blue", 20)),
     main="HSLike",ylab="theta")

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/HS-approx")
dev.copy2pdf(file = "hslike_profile.pdf")
dev.copy2eps(file = "hslike_profile.eps")

library(horseshoe)
par(mfrow=c(1,1))
data <-  seq(-5, 5, 0.05)
hsl_fit <- eval_hslike(data,nmc=5000)
tau.example <- HS.MMLE(data, 1)
#plot(data, HS.post.mean(data, tau = tau.example, Sigma2=mean(hsl_fit$Sigma2Save)),type="l",main="HS",ylab="theta")
plot(data, HS.post.mean(data, tau = tau.example, 1),type="l",main="HS",ylab="theta")
lines(data, colMeans(hsl_fit$ThetaSave), lty = 16,
     main="HSLike",ylab="theta")

