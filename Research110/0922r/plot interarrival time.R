alpha.est=apply(matrix(c(1,3,5),ncol=1),1,function(x) t3.log[5]+t3.log[6]*x)
sigma.est=apply(matrix(c(1,3,5),ncol=1),1,function(x) t3.nor[5]+t3.nor[6]*x)
beta.est=apply(matrix(c(1,3,5),ncol=1),1,function(x) t3.wei[5]+t3.wei[6]*x)
n = length(x[,1])


test.x <- seq(min(x), max(x), by = .0002)
cdf <- function(s){
  a <- s[1];b <- s[2];c <- s[3]
  plot(ecdf(as.numeric(x[,c(a,a+1,a+2)])),xlab = "Interarrival time",
       ylab = "Cumulative Proportion",main=paste("Distributions of interarrival time X under",c,"C"))
  lines(test.x,plnorm(test.x, meanlog=-alpha.est[b]^2/2,sdlog=alpha.est[b]),
        col="red")
  lines(test.x,pnorm(test.x, mean = 1, sd = sigma.est[b]),
        col="green")
  lines(test.x,pweibull(test.x,shape=beta.est[b],scale=(1/gamma(1/beta.est[b]+1))),
        col="blue")
  legend(0.985,1,                               
         col=1:4,  
         lwd=1,
         bty = "n",
         legend = c("data", 
                    paste("LogNormal (",paste(c(round(-alpha.est[b]^2/2,5),round(alpha.est[b],5)),collapse = ' ,'),")"),
                    paste("Normal (",paste(c(1,round(sigma.est[b],5)),collapse = ' , '),")"),
                    paste("Weibull (",paste(c(round((1/gamma(1/beta.est[b]+1)),5),round(beta.est[b],5)),collapse = ' , '),")")))
}


par(mfrow=c(1,3))
apply(matrix(c(1,4,1,2,1,3),ncol=3),1,cdf)
plot(ecdf(as.numeric(x[,7:8])),xlab = "Interarrival time",
     ylab = "Cumulative Proportion",main=paste("Distributions of interarrival time X under",5,"C"))
lines(test.x,plnorm(test.x, meanlog=-alpha.est[3]^2/2,sdlog=alpha.est[3]),
      col="red")
lines(test.x,pnorm(test.x, mean = 1, sd = sigma.est[3]),
      col="green")
lines(test.x,pweibull(test.x,shape=beta.est[3],scale=(1/gamma(1/beta.est[3]+1))),
      col="blue")
legend(0.96,1,                               
       col=1:4,  
       lwd=1,
       bty = "n",
       legend = c("data", 
                  paste("LogNormal (",paste(c(round(-alpha.est[3]^2/2,5),round(alpha.est[3],5)),collapse = ' ,'),")"),
                  paste("Normal (",paste(c(1,round(sigma.est[3],5)),collapse = ' , '),")"),
                  paste("Weibull (",paste(c(round((1/gamma(1/beta.est[3]+1)),5),round(beta.est[3],5)),collapse = ' , '),")")))

par(mfrow=c(1,1))

par(mfrow=c(3,3))
apply(matrix(1:8,ncol=1),1,function(s) plot(ecdf(x[,s])))
apply(matrix(1:8,ncol=1),1,function(s) plot((1:n - 1)/(n - 1), sort(x[,s]), type="l",
                                            main = print(paste(s,"th sample quantiles for the interarrival time of T-process")),
                                            xlab = "Sample Fraction",
                                            ylab = "Sample Quantile"))
par(mfrow=c(1,1))
library(EnvStats)
f <- function(s){
  a <- s[1];b <- s[2];c <- s[3]
  qqPlot(as.numeric(x[,c(a,a+1,a+2)]), dist = "lnorm", 
         param.list =list(meanlog=-alpha.est[b]^2/2,sdlog=alpha.est[b]), add.line = TRUE,
         xlab=paste("Quantiles of Normal (",paste(c('mean=','sd='),c(round(-alpha.est[b]^2/2,5),round(alpha.est[b],5)),collapse = ' , '),")"),
         ylab=paste("Quantiles of Log X under",c,"C"),
         main=paste("Normal Q-Q Plot for Log X under",c,"C"))
}
f.nor <- function(s){
  a <- s[1];b <- s[2];c <- s[3]
  qqPlot(as.numeric(x[,c(a,a+1,a+2)]), dist = "norm", 
         param.list =list(mean=1,sd=sigma.est[b]), add.line = TRUE,
         xlab=paste("Quantiles of Normal (",paste(c('mean=','sd='),c(1,round(sigma.est[b],5)),collapse = ', '),")"),
         ylab=paste("Quantiles of  X under",c,"C"),
         main=paste("Normal Q-Q Plot for  X under",c,"C"))
}
f.wei <- function(s){
  a <- s[1];b <- s[2];c <- s[3]
  qqPlot(as.numeric(x[,c(a,a+1,a+2)]), dist = "weibull", 
         param.list =list(shape=beta.est[b],scale=(1/gamma(1/beta.est[b]+1))), add.line = TRUE,
         xlab=paste("Quantiles of Weibull (",paste(c('scale=','shape='),c(round((1/gamma(1/beta.est[b]+1)),5),round(beta.est[b],5)),collapse = ', '),")"),
         ylab=paste("Quantiles of  X under",c,"C"),
         main=paste("Weibull Q-Q Plot for  X under",c,"C"))
}

par(mfrow=c(3,3))
apply(matrix(c(1,4,1,2,1,3),ncol=3),1,f)
qqPlot(as.numeric(x[,c(7,8)]), dist = "lnorm", 
       param.list =list(meanlog=-alpha.est[3]^2/2,sdlog=alpha.est[3]), add.line = TRUE,
       xlab=paste("Quantiles of Normal (",paste(c('mean=','sd='),c(round(-alpha.est[3]^2/2,5),round(alpha.est[3],5)),collapse = ', '),")"),
       ylab=paste("Quantiles of Log X under",5,"C"),
       main=paste("Normal Q-Q Plot for Log X under",5,"C"))
apply(matrix(c(1,4,1,2,1,3),ncol=3),1,f.nor)
qqPlot(as.numeric(x[,c(7,8)]), dist = "norm", 
       param.list =list(mean=1,sd=sigma.est[3]), add.line = TRUE,
       xlab=paste("Quantiles of Normal (",paste(c('mean=','sd='),c(1,round(sigma.est[3],5)),collapse = ', '),")"),
       ylab=paste("Quantiles of  X under",5,"C"),
       main=paste("Normal Q-Q Plot for  X under",5,"C"))
apply(matrix(c(1,4,1,2,1,3),ncol=3),1,f.wei)
qqPlot(as.numeric(x[,c(7,8)]), dist = "weibull", 
       param.list =list(shape=beta.est[3],scale=(1/gamma(1/beta.est[3]+1))), add.line = TRUE,
       xlab=paste("Quantiles of Weibull (",paste(c('scale=','shape='),c(round((1/gamma(1/beta.est[3]+1)),5),round(beta.est[3],5)),collapse = ' , '),")"),
       ylab=paste("Quantiles of  X under",5,"C"),
       main=paste("Weibull Q-Q Plot for  X under",5,"C"))
par(mfrow=c(1,1))






