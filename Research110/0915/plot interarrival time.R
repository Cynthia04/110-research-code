alpha.est=apply(matrix(c(1,3,5),ncol=1),1,function(x) t3.log[5]+t3.log[6]*x)
sigma.est=apply(matrix(c(1,3,5),ncol=1),1,function(x) t3.nor[5]+t3.nor[6]*x)
beta.est=apply(matrix(c(1,3,5),ncol=1),1,function(x) t3.wei[5]+t3.wei[6]*x)
n = length(x[,1])

par(mfrow=c(3,3))
apply(matrix(1:8,ncol=1),1,function(s) plot((1:n - 1)/(n - 1), sort(x[,s]), type="l",
                                            main = print(paste(s,"th sample quantiles for the interarrival time of transformed T-process")),
                                            xlab = "Sample Fraction",
                                            ylab = "Sample Quantile"))
par(mfrow=c(1,1))
f <- function(s){
  a <- s[1];b <- s[2]
  set.seed(11)
  t=cbind(sort(x[,a]),sort(rlnorm(45,meanlog=-alpha.est[b]^2/2,sdlog=alpha.est[b])),
          sort(rnorm(45,1,sigma.est[b])),
          sort(rweibull(45,shape=beta.est[b],scale=(1/gamma(1/beta.est[b]+1)))))
  matplot(x=(1:n - 1)/(n - 1),y=t,type="l",lty=1,lwd=1,col=1:4,xlab = "Sample Fraction",
          ylab = "Sample Quantile",main=print(paste(a,"th sample quantiles for the interarrival time of transformed T-process")))
  legend(0.7,1.0001,                                  
         col=1:4,  
         lwd=1,
         bty = "n",
         legend = c("data", "lognormal", "normal","weibull"))
  
}

par(mfrow=c(3,3))
apply(matrix(c(1:8,rep(1,3),rep(2,3),rep(3,2)),ncol=2),1,f)
par(mfrow=c(1,1))











