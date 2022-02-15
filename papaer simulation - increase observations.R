
##### papaer simulation #####

#---- 加大觀測個數

# table2 - F = N(1,sigma^2) 
est <- data.frame(a=c(9.219,9.269,9.235,9.422,9.503,9.434,9.649,9.638),
                  b=c(2.182,2.215,2.258,2.435,2.394,2.449,2.771,2.583),
                  sigma=c(2.339,2.080,1.782,3.357,2.518,2.541,3.944,4.053))

#n=100
n=100
set.seed(1109225017)
x <- apply(matrix(c(est[,3]*10^(-3)),ncol = 1),1,function(x){rnorm(n,mean=1,sd=x)})
x1 <- apply(x,2,function(s){cumsum(s)})
x2 <- matrix(0,ncol=3)
for(i in 1:8){
  t <- cbind (rep(est[i,1],n),rep(est[i,2],n),x1[,i])
  x2 <- rbind(x2,t)
}
x2 <- x2[-1,]
Ti<- apply(x2,1,function(x){log(x[2]*x[3]*10^(-2)/x[1]+1)/(x[2]*10^(-4))})
Ti <- matrix(Ti,ncol=8)
Ti <- rbind(rep(0,8),Ti)

#估參數
ll <- function(a,b,s,data){  #-loglikelihood
  sumt <- sum(data);J=length(data);sq <- s^2
  k <- c()
  for (s in 2:J){
    t<- (a/b*(exp(b*data[s])-exp(b*data[s-1]))-1)^2
    k <- append(k,t)
  }
  ans= -n/2*log(2*pi*sq)+n*log(a)+b*sumt-(sum(k)/(2*sq))
  -ans
}
ll(9.219*10^(-2),2.182*10^(-4),2.339*10^(-3),Ti[,1])
ll1 <- function(theta,data){
  ll(theta[1]*10^(-3),theta[2]*10^(-3),theta[3]*10^(-3),data)
}
ll1(c(9.219*10,2.182*10^(-1),2.339),Ti[,1])

set.seed(109225017);r1 <- runif(50,0,100);r2 <- runif(50,0,100);r3 <- runif(50,0,100)
apply(cbind(r1,r2,r3),1, FUN = function(x){ll1(x,Ti[,1])})
#apply(cbind(r1,r2,r3),1, FUN = function(x){optim(x,ll1,hessian=T,data=Ti[,1])})#確定無Inf or NA

table2_100 <- matrix(0,nrow = 8,ncol =7)
for(k in 1:8){
  optim.result=matrix(0,nrow = 1,ncol = 7)
  for (i in 1:50){
    fit=optim(c(r1[i],r2[i],r3[i]),ll1,hessian=T,data=Ti[,k])
    fit1=optim(fit$par,ll1,hessian=T,data=Ti[,k])
    if (abs(fit1$value-fit$value)>10^(-3)){
      fit <- fit1
      fit1<- optim(fit$par,ll1,hessian=T,data=Ti[,k])
    }
    I=solve(fit1$hessian)
    sda <- sqrt(I[1,1]*10^(-6))
    sdb <- sqrt(I[2,2]*10^(-6))
    sds <- sqrt(I[3,3]*10^(-6))
    t=cbind(fit1$par[1]*10^(-3),sda,
            fit1$par[2]*10^(-3),sdb,
            abs(fit1$par[3]*10^(-3)),sds,
            (-fit1$value))
    optim.result=rbind(optim.result,t)
  }
  optim.result <- optim.result[-1,]
  ind <- which.min(optim.result[,7])
  table2_100[k,] <- optim.result[ind,]
}

colnames(table2_100) <- c("a","sda","b","sdb","sigma","sds","logL")
table2_100

# table2_100 <- t(apply(table2_100,1,function(x){c(round(x[1]*10^2,3),x[2],
#                                          round(x[3]*10^4,3),x[4],
#                                          round(x[5]*10^3,3),x[6],
#                                          round(x[7],3))}))

#ci
ci_100 <- cbind(t(apply(table2_100[,1:2],1,function(x){x[1]+c(-1,1)*1.96*x[2]})),
                t(apply(table2_100[,3:4],1,function(x){x[1]+c(-1,1)*1.96*x[2]})),
                t(apply(table2_100[,5:6],1,function(x){x[1]+c(-1,1)*1.96*x[2]})))

colnames(ci_100) <- c("al","au","bl","bu","sl","su")
ci_100
?plotCI
#plot of CI
require(plotrix)
par(mfrow=c(1,3))
plotCI(1:8, est[,1]*10^(-2), li=ci_100[,1], ui=ci_100[,2],
       xlab="sample",ylab=expression(hat(a)),
       main="95% Confidence Interval of a")
plotCI(1:8, est[,2]*10^(-4), li=ci_100[,3], ui=ci_100[,4],
       ylim=c(0.00021,0.00027777),
       xlab="sample",ylab=expression(hat(b)),
       main="95% Confidence Interval of b")
plotCI(7,est[7,2]*10^(-4), li=ci_100[7,3], ui=ci_100[7,4])



plotCI(1:8, est[,3]*10^(-3), li=ci_100[,5], ui=ci_100[,6],
       xlab="sample",ylab=expression(hat(sigma)),
       main=expression(paste("95% Confidence Interval of ",sigma)))
par(mfrow=c(1,1))



#n=1000-----------------------------------------------------------------------------
n=1000
set.seed(1109225017)
x <- apply(matrix(c(est[,3]*10^(-3)),ncol = 1),1,function(x){rnorm(n,mean=1,sd=x)})
x1 <- apply(x,2,function(s){cumsum(s)})
x2 <- matrix(0,ncol=3)
for(i in 1:8){
  t <- cbind (rep(est[i,1],n),rep(est[i,2],n),x1[,i])
  x2 <- rbind(x2,t)
}
x2 <- x2[-1,]
Ti<- apply(x2,1,function(x){log(x[2]*x[3]*10^(-2)/x[1]+1)/(x[2]*10^(-4))})
Ti <- matrix(Ti,ncol=8)
Ti <- rbind(rep(0,8),Ti)

ll(9.219*10^(-2),2.182*10^(-4),2.339*10^(-3),Ti[,1])
ll1 <- function(theta,data){
  ll(theta[1]*10^(-3),theta[2]*10^(-3),theta[3]*10^(-3),data)
}
ll1(c(9.219*10,2.182*10^(-1),2.339),Ti[,1])

set.seed(109225017);r1 <- runif(50,0,100);r2 <- runif(50,0,1);r3 <- runif(50,0,10)
apply(cbind(r1,r2,r3),1, FUN = function(x){ll1(x,Ti[,1])})
apply(cbind(r1,r2,r3),1, FUN = function(x){optim(x,ll1,hessian=T,data=Ti[,1])})#確定無Inf or NA
try(optim(c(r1[1],r2[1],r3[1]),ll1,hessian=T,data=Ti[,1]))


table2_1000 <- matrix(0,nrow = 8,ncol =7)
for(k in 1:8){
  optim.result=matrix(0,nrow = 1,ncol = 7)
  for (i in 1:50){
    fit=try(optim(c(r1[i],r2[i],r3[i]),ll1,hessian=T,data=Ti[,k]))
    if(class(fit)=="try-error"){
      next
    }
    fit1=try(optim(fit$par,ll1,hessian=T,data=Ti[,k]))
    if(class(fit1)=="try-error"){
      next
    }
    
    if (abs(fit1$value-fit$value)>10^(-3)){
      fit <- fit1
      fit1<- optim(fit$par,ll1,hessian=T,data=Ti[,k])
    }
    I=solve(fit1$hessian)
    sda <- sqrt(I[1,1]*10^(-6))
    sdb <- sqrt(I[2,2]*10^(-6))
    sds <- sqrt(I[3,3]*10^(-6))
    t=cbind(fit1$par[1]*10^(-3),sda,
            fit1$par[2]*10^(-3),sdb,
            abs(fit1$par[3]*10^(-3)),sds,
            (-fit1$value))
    optim.result=rbind(optim.result,t)
  }
  optim.result <- optim.result[-1,]
  ind <- which.min(optim.result[,7])
  table2_1000[k,] <- optim.result[ind,]
  
}

colnames(table2_1000) <- c("a","sda","b","sdb","sigma","sds","logL")
table2_1000

# table2_1000 <- t(apply(table2_1000,1,function(x){c(round(x[1]*10^2,3),x[2],
#                                          round(x[3]*10^4,3),x[4],
#                                          round(x[5]*10^3,3),x[6],
#                                          round(x[7],3))}))

#ci
ci_1000 <- cbind(t(apply(table2_1000[,1:2],1,function(x){x[1]+c(-1,1)*1.96*x[2]})),
                 t(apply(table2_1000[,3:4],1,function(x){x[1]+c(-1,1)*1.96*x[2]})),
                 t(apply(table2_1000[,5:6],1,function(x){x[1]+c(-1,1)*1.96*x[2]})))

colnames(ci_1000) <- c("al","au","bl","bu","sl","su")
ci_1000

#plot of CI
require(plotrix)
par(mfrow=c(1,3))
plotCI(1:8, est[,1]*10^(-2), li=ci_1000[,1], ui=ci_1000[,2],
       xlab="sample",ylab=expression(hat(a)),
       main="95% Confidence Interval of a")
plotCI(1:8, est[,2]*10^(-4), li=ci_1000[,3], ui=ci_1000[,4],
       xlab="sample",ylab=expression(hat(b)),
       main="95% Confidence Interval of b")
plotCI(1:8, est[,3]*10^(-3), li=ci_1000[,5], ui=ci_1000[,6],
       xlab="sample",ylab=expression(hat(sigma)),
       main=expression(paste("95% Confidence Interval of ",sigma)))
par(mfrow=c(1,1))


