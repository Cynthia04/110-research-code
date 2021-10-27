#------
#0831
#Q4
#------
#normal
t3 <- matrix(c(9.15 *10^(-2),
               1.01 *10^(-3),
               2.07 *10^(-4),
               1.19 *10^(-5),
               2.76 *10^(-3),
               4.24 *10^(-4)),ncol=2,byrow = T)

est <- apply(t3,1,function(x) x[1]+c(1,3,5)*x[2])
colnames(est) <- c("a","b","s")
rownames(est) <- c("1C","3C","5C")
est0 <- rbind(matrix(rep(est[1,],3),ncol=3,byrow = T),
              matrix(rep(est[2,],3),ncol=3,byrow = T),
              matrix(rep(est[3,],2),ncol=3,byrow = T))

colnames(est0) <- c("a","b","s")


set.seed(1109225017)
x <- apply(matrix(est0[,3],ncol=1),1,
           function(x){rnorm(45,mean = 1,sd=x)})
x1 <- apply(x,2,function(s){cumsum(s)})
x2 <- matrix(0,ncol=3)
for(i in 1:8){
  t <- cbind (rep(est0[i,1],45),rep(est0[i,2],45),x1[,i])
  x2 <- rbind(x2,t)
}
x2 <- x2[-1,]
Ti<- apply(x2,1,function(x){log(x[2]*x[3]/x[1]+1)/x[2]})
Ti <- matrix(Ti,ncol=8)
Ti <- rbind(rep(0,8),Ti)
ll <- function(a,b,s,data){  #-loglikelihood
  sumt <- sum(data);J=length(data);sq <- s^2
  k <- c()
  for (s in 2:J){
    t<- (a/b*(exp(b*data[s])-exp(b*data[s-1]))-1)^2
    k <- append(k,t)
  }
  ans= -45/2*log(2*pi*sq)+45*log(a)+b*sumt-(sum(k)/(2*sq))
  -ans
}

ll1 <- function(theta,data){
  ll(theta[1]*10^(-3),theta[2]*10^(-3),theta[3]*10^(-3),data)
}

set.seed(109225017);r1 <- runif(50,0,100);r2 <- runif(50,0,100);r3 <- runif(50,0,100)
table2 <- matrix(0,nrow = 8,ncol =4)
system.time({
  for(k in 1:8){
    optim.result=matrix(0,nrow = 1,ncol = 4)
    for (i in 1:50){
      fit=optim(c(r1[i],r2[i],r3[i]),ll1,hessian=T,data=Ti[,k])
      fit1=optim(fit$par,ll1,hessian=T,data=Ti[,k])
      while (abs(fit1$value-fit$value)>10^(-3)){
        fit <- fit1
        fit1<- optim(fit$par,ll1,hessian=T,data=Ti[,k])
      }
      t=cbind(fit1$par[1]*10^(-3),
              fit1$par[2]*10^(-3),
              abs(fit1$par[3]*10^(-3)),
              (-fit1$value))
      optim.result=rbind(optim.result,t)
    }
    optim.result <- optim.result[-1,]
    print(optim.result)
    ind <- which.max(optim.result[,4])
    table2[k,] <- optim.result[ind,]
  }
})
colnames(table2) <- c("a","b","sigma","logL")
table2.nor <- table2

ll.log <- function(a,b,c,data){  #-loglikelihood
  sumt <- sum(data);J=length(data);
  k <- c()
  for (s in 2:J){
    t<- log(exp(b*data[s])-exp(b*data[s-1]))
    k <- append(k,t)
  }
  k1 <- c()
  for (s in 2:J){
    t1<- (log(a)+log(exp(b*data[s])-exp(b*data[s-1]))-log(b)+c^2/2)^2
    k1 <- append(k1,t1)
  }
  ans= 45*log(b)-sum(k)-45*log(c)-45/2*log(2*pi)-1/(2*c^2)*sum(k1)+b*sumt
  -ans
}
ll1.log <- function(theta,data){
  ll.log(theta[1]*10^(-2),theta[2]*10^(-4),theta[3]*10^(-3),data)
}
set.seed(1092501);r1 <- runif(50,0,10);r2 <- runif(50,0,10);r3 <- runif(50,0,10)
table2 <- matrix(0,nrow = 8,ncol =4)
system.time({
  for(k in 1:8){
    optim.result=matrix(0,nrow = 1,ncol = 4)
    for (i in 1:50){
      fit=optim(c(r1[i],r2[i],r3[i]),ll1.log,hessian=T,data=Ti[,k])
      fit1=optim(fit$par,ll1.log,hessian=T,data=Ti[,k])
      while (abs(fit1$value-fit$value)>10^(-3)){
        fit <- fit1
        fit1<- optim(fit$par,ll1.log,hessian=T,data=Ti[,k])
      }
      t=cbind(fit1$par[1]*10^(-2),
              fit1$par[2]*10^(-4),
              abs(fit1$par[3]*10^(-3)),
              (-fit1$value))
      optim.result=rbind(optim.result,t)
    }
    optim.result <- optim.result[-1,]
    print(optim.result)
    ind <- which.max(optim.result[,4])
    table2[k,] <- optim.result[ind,]
  }
})
colnames(table2) <- c("a","b","alpha","logL")
table2.log <- table2

ll.wei = function(a,b,c,data){  
  sumt =sum(data);J=length(data)
  k <- c()
  for (s in 2:J){
    t<-  (c-1)*(lgamma(1/c+1)+log((a/b)*(exp(b*data[s])-exp(b*data[s-1]))))
    k <- append(k,t)
  }
  k1 <- c()
  for (s in 2:J){
    t1<- (a/b*gamma(1/c+1)*(exp(b*data[s])-exp(b*data[s-1])))^c
    k1 <- append(k1,t1)
  }
  ans= 45*lgamma(1/c+1)+45*log(c)+45*log(a)+sum(k)-sum(k1)+b*sumt
  -ans
}
ll1.wei <- function(theta,data){
  ll.wei(theta[1]*10^(-2),theta[2]*10^(-4),theta[3]*10^2,data)
}
set.seed(10922501);r1 <- runif(50,0,10);r2 <- runif(50,0,10);r3 <- runif(50,0,10)
apply(cbind(r1,r2,r3),1, FUN = function(x){ll1.wei(x,Ti[,1])})
table2 <- matrix(0,nrow = 8,ncol =4)
system.time({
  for(k in 1:8){
    optim.result=matrix(0,nrow = 1,ncol = 4)
    for (i in 1:50){
      fit=optim(c(r1[i],r2[i],r3[i]),ll1.wei,hessian=T,data=Ti[,k])
      fit1=optim(fit$par,ll1.wei,hessian=T,data=Ti[,k])
      while (abs(fit1$value-fit$value)>10^(-3)){
        fit <- fit1
        fit1<- optim(fit$par,ll1.wei,hessian=T,data=Ti[,k])
      }
      t=cbind(fit1$par[1]*10^(-2),
              fit1$par[2]*10^(-4),
              fit1$par[3]*10^(2),
              (-fit1$value))
      optim.result=rbind(optim.result,t)
    }
    optim.result <- optim.result[-1,]
    print(optim.result)
    ind <- which.max(optim.result[,4])
    table2[k,] <- optim.result[ind,]
  }
})
colnames(table2) <- c("a","b","beta","logL")
table2.wei <- table2
result <- list(log=table2.log,
               wei=table2.wei,
               nor=table2.nor)
library(data.table)
outputfile <- "result.csv" #output file name
sep <- "," #define the separator (related to format of the output file)
for(i in names(result)){
  fwrite(list(i), file=outputfile, sep=sep, append=T) #write names of the list elements
  ele <- result[[i]]
  fwrite(ele, file=outputfile, sep=sep, append=T, col.names=T) 
  fwrite(list(NA), file=outputfile, append=T) #add an empty row to separate elements
}
#------------------------------------------------------------------------------------
#table3
#------------------------------------------------------------------------------------

tll <- function(a0,a1,b0,b1,s0,s1,data){ # ll is -loglikeihood, so here is also -logL
  aS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) a0+a1*x)
  bS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) b0+b1*x)
  sS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) s0+s1*x)
  t1 <- sum(apply(matrix(1:3,ncol = 1),1,function(x) ll(aS[1],bS[1],sS[1],data[,x])))
  t2 <- sum(apply(matrix(4:6,ncol = 1),1,function(x) ll(aS[2],bS[2],sS[2],data[,x])))
  t3 <- sum(apply(matrix(7:8,ncol = 1),1,function(x) ll(aS[3],bS[3],sS[3],data[,x])))
  ans <- t1+t2+t3
  ans
}
tll1 <- function(theta,data){
  tll(theta[1]*10^(-2),theta[2]*10^(-3),theta[3]*10^(-4),theta[4]*10^(-5),theta[5]*10^(-3),theta[6]*10^(-4),data)
}
set.seed(1);r1 <- runif(50,0,100);r2 <- runif(50,0,100);r3 <- runif(50,0,100);r4 <- runif(50,0,100);r5 <- runif(50,0,100);r6 <- runif(50,0,100)
optim.result=matrix(0,nrow = 1,ncol = 7)
for (i in 1:50){
  fit=try(optim(c(r1[i],r2[i],r3[i],r4[i],r5[i],r6[i]),tll1,hessian=T,data=Ti))
  if(class(fit)=="try-error"){
    next
  }
  fit1=try(optim(fit$par,tll1,hessian=T,data=Ti))
  if(class(fit1)=="try-error"){
    next
  }
  while (abs(fit1$value-fit$value)>10^(-3)){
    fit <- fit1
    fit1<- optim(fit$par,tll1,hessian=T,data=Ti)
  }
  t=cbind(fit1$par[1]*10^(-2),fit1$par[2]*10^(-3),fit1$par[3]*10^(-4),
          fit1$par[4]*10^(-5),abs(fit1$par[5]*10^(-3)),abs(fit1$par[6]*10^(-4)),
          (-fit1$value))
  optim.result=rbind(optim.result,t)
}
optim.result <- optim.result[-1,]
ind <- which.max(optim.result[,7])
t3.nor <- matrix(optim.result[ind,],ncol = 1)

t3sim <- t3.nor[1:6]
S0 <- c()
for (i in c(1,3,5)){
  t <- t3sim[i]+0.5*t3sim[i+1]
  S0 <- append(S0,t)
}
S0
EZ <- function(i){
  a <- S0[1];b <- S0[2];s <- S0[3]
  k <- a/b
  1/b*(log((i+k)/(i-1+k))+s^2/2*((i-1)/(i-1+k)^2-i/(i+k)^2))
}
EZ(1)

x=1;ans <- EZ(x)
while( ans>8){
  x=x+1 
  ans=EZ(x)
}
x
ans



tll.log <- function(a0,a1,b0,b1,s0,s1,data){ # ll is -loglikeihood, so here is also -logL
  aS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) a0+a1*x)
  bS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) b0+b1*x)
  sS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) s0+s1*x)
  t1 <- sum(apply(matrix(1:3,ncol = 1),1,function(x) ll.log(aS[1],bS[1],sS[1],data[,x])))
  t2 <- sum(apply(matrix(4:6,ncol = 1),1,function(x) ll.log(aS[2],bS[2],sS[2],data[,x])))
  t3 <- sum(apply(matrix(7:8,ncol = 1),1,function(x) ll.log(aS[3],bS[3],sS[3],data[,x])))
  ans <- t1+t2+t3
  ans
}
tll1.log <- function(theta,data){
  tll.log(theta[1]*10^(-2),theta[2]*10^(-3),theta[3]*10^(-4),theta[4]*10^(-5),theta[5]*10^(-3),theta[6]*10^(-4),data)
}
set.seed(1);r1 <- runif(50,0,100);r2 <- runif(50,0,100);r3 <- runif(50,0,100);r4 <- runif(50,0,100);r5 <- runif(50,0,100);r6 <- runif(50,0,100)
optim.result=matrix(0,nrow = 1,ncol = 7)
for (i in 1:50){
  fit=try(optim(c(r1[i],r2[i],r3[i],r4[i],r5[i],r6[i]),tll1.log,hessian=T,data=Ti))
  if(class(fit)=="try-error"){
    next
  }
  fit1=try(optim(fit$par,tll1.log,hessian=T,data=Ti))
  if(class(fit1)=="try-error"){
    next
  }
  while (abs(fit1$value-fit$value)>10^(-3)){
    fit <- fit1
    fit1<- optim(fit$par,tll1.log,hessian=T,data=Ti)
  }
  t=cbind(fit1$par[1]*10^(-2),fit1$par[2]*10^(-3),fit1$par[3]*10^(-4),
          fit1$par[4]*10^(-5),abs(fit1$par[5]*10^(-3)),abs(fit1$par[6]*10^(-4)),
          (-fit1$value))
  optim.result=rbind(optim.result,t)
}
optim.result <- optim.result[-1,]
ind <- which.max(optim.result[,7])
t3.log <- matrix(optim.result[ind,],ncol=1)

t3sim <- t3.log[1:6]
S0 <- c()
for (i in c(1,3,5)){
  t <- t3sim[i]+0.5*t3sim[i+1]
  S0 <- append(S0,t)
}
S0
EZ <- function(i){
  a <- S0[1];b <- S0[2];c <- S0[3]
  k <- a/b
  1/b*(log((i+k)/(i-1+k))+(exp(1)^(c^2)-1)/2*((i-1)/(i-1+k)^2-i/(i+k)^2))
}
EZ(1)

x=1;ans <- EZ(x)
while( ans>8){
  x=x+1 
  ans=EZ(x)
}
x
ans


tll.wei<- function(a0,a1,b0,b1,s0,s1,data){ # ll is -loglikeihood, so here is also -logL
  aS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) a0+a1*x)
  bS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) b0+b1*x)
  sS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) s0+s1*x)
  t1 <- sum(apply(matrix(1:3,ncol = 1),1,function(x) ll.wei(aS[1],bS[1],sS[1],data[,x])))
  t2 <- sum(apply(matrix(4:6,ncol = 1),1,function(x) ll.wei(aS[2],bS[2],sS[2],data[,x])))
  t3 <- sum(apply(matrix(7:8,ncol = 1),1,function(x) ll.wei(aS[3],bS[3],sS[3],data[,x])))
  ans <- t1+t2+t3
  ans
}
tll1.wei <- function(theta,data){
  tll.wei(theta[1]*10^(-2),theta[2]*10^(-3),theta[3]*10^(-4),theta[4]*10^(-5),theta[5],theta[6],data)
}
N=50
rn <- matrix(0,ncol=7,nrow=N)
set.seed(1);rn[,1]<-runif(N,0,10);rn[,2]<-runif(N,0,10);rn[,3]<-runif(N,0,10)
set.seed(1);rn[,4]<-runif(N,0,10);rn[,5]<-runif(N,0,10);rn[,6]<-runif(N,0,10)
rn[,7]<-1:N
Opt<-function(x){
  iter=x[7]
  initial=c(x[1],x[2],x[3],x[4],x[5],x[6])
  error=1
  before=0
  while(abs(error)-1e-15>0){
    step1=try(optim(initial,tll1.wei,control=list(factr=1e-15,maxit=10^9),data=Ti),silent = T)
    if(step1$convergence==0){
      initial <- step1$par
      error <- step1$value-before
      before <- step1$value
    }
  }
  return(c(step1$par,-step1$value)) 
}
result=apply(as.matrix(rn),1,function(x){Opt(x)})
result <- t(result)
ind1 <- which.max(result[,7])
t3.wei <- matrix(c(result[ind1,1]*10^(-2),result[ind1,2]*10^(-3),result[ind1,3]*10^(-4),
                   result[ind1,4]*10^(-5),result[ind1,5],result[ind1,6],result[ind1,7]),ncol=1)

t3sim <- t3.wei[1:6]
S0 <- c()
for (i in c(1,3,5)){
  t <- t3sim[i]+0.5*t3sim[i+1]
  S0 <- append(S0,t)
}
S0
EZ <- function(i){
  a <- S0[1];b <- S0[2];c <- S0[3]
  k <- a/b; var <- (1/gamma(1/c+1))^2*gamma(2/c+1)-1
  1/b*(log((i+k)/(i-1+k))+var/2*((i-1)/(i-1+k)^2-i/(i+k)^2))
}
EZ(1)

x=1;ans <- EZ(x)
while( ans>8){
  x=x+1 
  ans=EZ(x)
}
x
ans



result <- list(Lognormal=t3.log,
               Weibull=t3.wei,
               Normal=t3.nor)
library(data.table)
outputfile <- "result.csv" #output file name
sep <- "," #define the separator (related to format of the output file)
for(i in names(result)){
  fwrite(list(i), file=outputfile, sep=sep, append=T) #write names of the list elements
  ele <- result[[i]]
  fwrite(ele, file=outputfile, sep=sep, append=T, col.names=T) 
  fwrite(list(NA), file=outputfile, append=T) #add an empty row to separate elements
}



#-------------------------------------------------------------------------------
# 5.1 Model checking
#-------------------------------------------------------------------------------
#normal
#residuals
t3sim <- t3.log[1:6,1]
est.re <- data.frame(a=apply(matrix(c(1,3,5),ncol = 1),1,function(x) t3sim[1]+t3sim[2]*x),
                     b=apply(matrix(c(1,3,5),ncol = 1),1,function(x) t3sim[3]+t3sim[4]*x),
                     sigma=apply(matrix(c(1,3,5),ncol = 1),1,function(x) t3sim[5]+t3sim[6]*x))
#----------------- Simulate data m_jk=45 (n_1,n_2,n_3)=(3,3,2) ----------------#
set.seed(1109225017)
x.re <- apply(matrix(c(est.re[,3]),ncol = 1),1,function(x){rlnorm(45,meanlog=-x^2/2,sdlog=x^2)})
x1.re <- apply(x.re,2,function(s){cumsum(s)})
x2.re <- matrix(0,ncol=3)
for(i in 1:3){
  t <- cbind (rep(est.re[i,1],45),rep(est.re[i,2],45),x1.re[,i])
  x2.re <- rbind(x2.re,t)
}
x2.re <- x2.re[-1,]
Ti.re<- apply(x2.re,1,function(x){log(x[2]*x[3]/x[1]+1)/(x[2])})
Ti.re <- matrix(Ti.re,ncol=3)
Ti.re <- rbind(rep(0,3),Ti.re)


simulationdata1 <- Ti[2:46,]-Ti[1:45,]
matplot(apply(matrix(1:3,ncol = 1),1,function(x) Ti[,x]-Ti.re[,1]))

simulationdata1.re <- Ti.re[2:46,]-Ti.re[1:45,]

par(mfrow=c(1,3))
matplot(apply(matrix(1:3,ncol = 1),1,function(x) simulationdata1[,x]-simulationdata1.re[,1]),
        xlab = "g-cycle",ylab="Residuals",main="Residual plot in 1C",
        pch=1,col="black")
abline(h=0)

matplot(apply(matrix(4:6,ncol = 1),1,function(x) simulationdata1[,x]-simulationdata1.re[,2]),
        xlab = "g-cycle",ylab="Residuals",main="Residual plot in 3C",
        pch=1,col="black")
abline(h=0)

matplot(apply(matrix(7:8,ncol = 1),1,function(x) simulationdata1[,x]-simulationdata1.re[,3]),
        xlab = "g-cycle",ylab="Residuals",main="Residual plot in 5C",
        pch=1,col="black")
abline(h=0)
par(mfrow=c(1,1))


#normal test

nt <- cbind(est.re[,1]/est.re[,2],est.re[,2])

sq1 <- apply(matrix(1:3,ncol=1),1,function(x) nt[1,1]*(exp(nt[1,2]*Ti[2:46,x])-exp(nt[1,2]*Ti[1:45,x])))
sq3 <- apply(matrix(4:6,ncol=1),1,function(x) nt[2,1]*(exp(nt[2,2]*Ti[2:46,x])-exp(nt[2,2]*Ti[1:45,x])))
sq5 <- apply(matrix(7:8,ncol=1),1,function(x) nt[3,1]*(exp(nt[3,2]*Ti[2:46,x])-exp(nt[3,2]*Ti[1:45,x])))

par(mfrow=c(1,3))
qqnorm(as.numeric(sq1), pch = 1, frame = FALSE)
qqline(as.numeric(sq1), col = "steelblue", lwd = 2)

qqnorm(as.numeric(sq3), pch = 1, frame = FALSE)
qqline(as.numeric(sq3), col = "steelblue", lwd = 2)

qqnorm(as.numeric(sq5), pch = 1, frame = FALSE)
qqline(as.numeric(sq5), col = "steelblue", lwd = 2)
par(mfrow=c(1,1))


ks.test(as.numeric(sq1),rnorm(135,1,est.re[1,3]))
ks.test(as.numeric(sq3),rnorm(135,1,est.re[2,3]))
ks.test(as.numeric(sq5),rnorm(90,1,est.re[3,3]))

#independent 
sq <- cbind(sq1,sq3,sq5)
library(Kendall)
apply(matrix(1:8,ncol=1),1,function(x) MannKendall(sq[,x]))

#lognormal


