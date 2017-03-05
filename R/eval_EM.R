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