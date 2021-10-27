#--------------- Use the normal estimate of table 2 to simulate data -----------#
est <- data.frame(a=c(9.219,9.269,9.235,9.422,9.503,9.434,9.649,9.638),
                  b=c(2.182,2.215,2.258,2.435,2.394,2.449,2.771,2.583),
                  sigma=c(2.339,2.080,1.782,3.357,2.518,2.541,3.944,4.053))
#----------------- Simulate data m_jk=45 (n_1,n_2,n_3)=(3,3,2) ----------------#
set.seed(1109225017)
x <- apply(matrix(c(est[,3]*10^(-3)),ncol = 1),1,function(x){rnorm(45,mean=1,sd=x)})
x1 <- apply(x,2,function(s){cumsum(s)})
x2 <- matrix(0,ncol=3)
for(i in 1:8){
  t <- cbind (rep(est[i,1],45),rep(est[i,2],45),x1[,i])
  x2 <- rbind(x2,t)
}
x2 <- x2[-1,]
Ti<- apply(x2,1,function(x){log(x[2]*x[3]*10^(-2)/x[1]+1)/(x[2]*10^(-4))})
Ti <- matrix(Ti,ncol=8)
Ti <- rbind(rep(0,8),Ti)
rm(i,t)

matplot(Ti,type="b",pch=1:8,col=1:8,
        xlab = "g-cycle",ylab = expression(CR[g]),
        main="Recurrent event proceses for the battery simulation data")

#----------------------------- plot Figure 2 ------------------------------#

simulationdata1 <- Ti[2:46,]-Ti[1:45,]
asd <- simulationdata1

simud1.1 <- data.frame(y=as.vector(asd[,1:3]),
                       x=rep(1:45,3))
simud1.3 <- data.frame(y=as.vector(asd[,4:6]),
                       x=rep(1:45,3))
simud1.5 <- data.frame(y=as.vector(asd[,7:8]),
                       x=rep(1:45,2))
lm1.1 <- lm(y~x,data=simud1.1)
lm1.3 <- lm(y~x,data=simud1.3)
lm1.5 <- lm(y~x,data=simud1.5)

matplot(asd,xlim = c(0,45),ylim = c(9.0,11.0),
        xlab = "g-cycle",ylab = expression(CR[g]),
        main=expression(paste("Interarrival times of events ",list(CR[g]),"s")),
        pch=c(1,1,1,2,2,2,4,4),col="black")

abline(lm1.1,lwd=3,col="green")
abline(lm1.3,lwd=3,col="#63B8FF",lty=5)
abline(lm1.5,lwd=3,col="red1",lty=4)

legend(38,11,                               
       pch = c(1,2,4),  
       bty = "n",
       cex = 0.7,
       legend = c("1C observation", "3C observation", "5C observation"))

legend(0,9.6,                               
       col = c("green","#63B8FF","red1"),  
       lwd=3,
       lty=c(1,5,4),
       bty = "n",
       cex = 0.7,
       legend = c("1C fitted line", "3C fitted line", "5C fitted line"))

#--------------------------- simulate table2 ----------------------------------#

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

table2 <- matrix(0,nrow = 8,ncol =7)
system.time({
  for(k in 1:8){
    optim.result=matrix(0,nrow = 1,ncol = 7)
    for (i in 1:50){
      fit=optim(c(r1[i],r2[i],r3[i]),ll1,hessian=T,data=Ti[,k])
      fit1=optim(fit$par,ll1,hessian=T,data=Ti[,k])
      while (abs(fit1$value-fit$value)>10^(-3)){
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
    print(optim.result)
    ind <- which.max(optim.result[,7])
    table2[k,] <- optim.result[ind,]
  }
})
colnames(table2) <- c("a","sda","b","sdb","sigma","sds","logL")
# table2 <- t(apply(table2,1,function(x){c(round(x[1]*10^2,3),x[2],
#                              round(x[3]*10^4,3),x[4],
#                              round(x[5]*10^3,3),x[6],
#                              round(x[7],3))}))
# 
table2.nor <- table2[,c(1,3,5,7)]



#Figure 5

S.abs <- data.frame(S=table2[,c(1,3,5)],
                    group=c(1,1,1,3,3,3,5,5))
lma <- lm(S.a~group,data=S.abs[,c(1,4)])
lmb <- lm(S.b~group,data=S.abs[,c(2,4)])
lms <- lm(S.sigma~group,data=S.abs[,c(3,4)])


par(mfrow=c(1,3))
plot(S.abs$group,S.abs$S.a,xlab = "Discharge Rate",ylab = "a",main="normal")
abline(lma)

plot(S.abs$group,S.abs$S.b,xlab = "Discharge Rate",ylab = "b",main="normal")
abline(lmb)

plot(S.abs$group,S.abs$S.sigma,xlab = "Discharge Rate",
     ylab = expression(sigma),main="normal")
abline(lms)
par(mfrow=c(1,1))

table3 <- rbind(summary(lma)$coef[,1:2],summary(lmb)$coef[,1:2],summary(lms)$coef[,1:2])
row.names(table3) <- c(expression(a_0),expression(a_1),
                       expression(b_0),expression(b_1),
                       expression(sigma_0),expression(sigma_1))

#------------------- use regression method to  sumulate table3 ------------------#

table3.nor <- cbind(table3,rbind(confint(lma),confint(lmb),confint(lms)))

#------------------- use optim to  sumulate table3 -------------------------------#
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
apply(cbind(r1,r2,r3,r4,r5,r6),1, FUN = function(x){tll1(x,data=Ti)})
#apply(cbind(r1,r2,r3,r4,r5,r6),1, FUN = function(x){optim(x,tll1,hessian=F,data=Ti)})
start_time <- Sys.time()
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
end_time <- Sys.time()
end_time - start_time
optim.result <- optim.result[-1,]
ind <- which.max(optim.result[,7])
t3.nor <- matrix(optim.result[ind,],ncol = 1)
t3sim <- optim.result[ind,1:6]
tll( 0.09146606,0.001014137,0.0002073602, 1.189044e-05,0.002759585,0.000423787,Ti)
tll( 0.091488761,0.001004613,0.000207317,1.20E-05,0.00172493,0.000346193,Ti)



#---------------------
#          MLE
#---------------------
# a0   0.09146606
# a1   0.001014137
# b0   0.0002073602
# b1   1.189044e-05
# s0   0.002759585
# s1   0.000423787
#logL  6.588330e+02
#---------------------
t3sim[1]*10^2;t3sim[2]*10^3;t3sim[3]*10^4;t3sim[4]*10^5;t3sim[5]*10^3;t3sim[6]*10^4
#---------------------
#      MLE(同論文 -> 四捨五入到小數點第二位)
#---------------------
# a0   9.15*10^(-2)
# a1   1.01*10^(-3)
# b0   2.07*10^(-4)
# b1   1.19*10^(-5)
# s0   2.76*10^(-3)
# s1   4.24*10^(-4)
#logL  6.588330e+02
#---------------------

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
#-------------------------------------------------------------------------------
# 5.1 Model checking
#-------------------------------------------------------------------------------
#residuals

Ti[,1:3]

simulationdata1[,1:3]
t3sim


est.re <- data.frame(a=apply(matrix(c(1,3,5),ncol = 1),1,function(x) t3sim[1]+t3sim[2]*x),
                  b=apply(matrix(c(1,3,5),ncol = 1),1,function(x) t3sim[3]+t3sim[4]*x),
                  sigma=apply(matrix(c(1,3,5),ncol = 1),1,function(x) t3sim[5]+t3sim[6]*x))
#----------------- Simulate data m_jk=45 (n_1,n_2,n_3)=(3,3,2) ----------------#
set.seed(1109225017)
x.re <- apply(matrix(c(est.re[,3]),ncol = 1),1,function(x){rnorm(45,mean=1,sd=x)})
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


#-------------------------------------------------------------------------------
# 6 Compare of EOP
#-------------------------------------------------------------------------------

plot(1:45,simulationdata1[,1])







#-----------------------------------------------------------------------------
#lognormal 
#-----------------------------------------------------------------------------
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

ll.log(9.219*10^(-2),2.182*10^(-4),2.34*10^(-3),Ti[,1])


ll1.log <- function(theta,data){
  ll.log(theta[1]*10^(-3),theta[2]*10^(-3),theta[3]*10^(-3),data)
}
set.seed(10922501);r1 <- runif(50,0,100);r2 <- runif(50,0,100);r3 <- runif(50,0,100)
table2 <- matrix(0,nrow = 8,ncol =7)
system.time({
  for(k in 1:8){
    optim.result=matrix(0,nrow = 1,ncol = 7)
    for (i in 1:50){
      fit=optim(c(r1[i],r2[i],r3[i]),ll1.log,hessian=T,data=Ti[,k])
      fit1=optim(fit$par,ll1.log,hessian=T,data=Ti[,k])
      while (abs(fit1$value-fit$value)>10^(-3)){
        fit <- fit1
        fit1<- optim(fit$par,ll1.log,hessian=T,data=Ti[,k])
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
    print(optim.result)
    ind <- which.max(optim.result[,7])
    table2[k,] <- optim.result[ind,]
  }
})
colnames(table2) <- c("a","sda","b","sdb","sigma","sds","logL")
table2.log <- table2[,c(1,3,5,7)]

#-------------------------------------------------------------------------------
#Figure 5
#-------------------------------------------------------------------------------

S.log <- data.frame(S=table2[,c(1,3,5)],
                    group=c(1,1,1,3,3,3,5,5))
lma.log <- lm(S.a~group,data=S.log[,c(1,4)])
lmb.log <- lm(S.b~group,data=S.log[,c(2,4)])
lms.log <- lm(S.sigma~group,data=S.log[,c(3,4)])


par(mfrow=c(1,3))
plot(S.log$group,S.log$S.a,xlab = "Discharge Rate",ylab = "a",main="log-normal")
abline(lma.log)

plot(S.log$group,S.log$S.b,xlab = "Discharge Rate",ylab = "b",main="log-normal")
abline(lmb.log)

plot(S.log$group,S.log$S.sigma,xlab = "Discharge Rate",
     ylab = expression(alpha),main="log-normal")
abline(lms.log)
par(mfrow=c(1,1))
table3.log <- rbind(summary(lma.log)$coef[,1:2],summary(lmb.log)$coef[,1:2],summary(lms.log)$coef[,1:2])
row.names(table3.log) <- c(expression(a_0),expression(a_1),
                       expression(b_0),expression(b_1),
                       expression(alpha_0),expression(alpha_1))

#------------------- use regression method to  sumulate table3 ------------------#

table3.log <- cbind(table3.log,rbind(confint(lma.log),confint(lmb.log),confint(lms.log)))

#-------------------------------------------------------------------------------

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
apply(cbind(r1,r2,r3,r4,r5,r6),1, FUN = function(x){tll1.log(x,data=Ti)})
#apply(cbind(r1,r2,r3,r4,r5,r6),1, FUN = function(x){optim(x,tll1,hessian=F,data=Ti)})
start_time <- Sys.time()
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
end_time <- Sys.time()
end_time - start_time
optim.result <- optim.result[-1,]
ind <- which.max(optim.result[,7])
t3.log <- matrix(optim.result[ind,],ncol=1)
tll.log(0.091489121,
        0.001004538,
        0.000207301,
        1.20E-05,
        0.001724255,
        0.000346295,Ti)

#table3
#---------------------
#          MLE
#---------------------
# a0   0.09146624
# a1   0.001013941
# b0   0.0002073534
# b1   1.189871e-05
# s0   0.002758018
# s1   0.0004241077
#logL  6.589139e+02
#---------------------
#---------------------
#      MLE(同論文 -> 四捨五入到小數點第二位)
#---------------------
# a0   9.15*10^(-2)
# a1   1.01*10^(-3)
# b0   2.07*10^(-4)
# b1   1.19*10^(-5)
# s0   2.76*10^(-3)
# s1   4.24*10^(-4)
#logL  6.589139e+02
#---------------------





#-----------------------------------------------------------------------------
#Weibull 


#-----------------------------------------------------------------------------
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
ll.wei(9.222*10^(-2),2.164*10^(-4),5.565*10^(2),Ti[,1])
ll1.wei <- function(theta,data){
  ll.wei(theta[1]*10^(-2),theta[2]*10^(-4),theta[3]*10^2,data)
}
set.seed(10922501);r1 <- runif(50,0,10);r2 <- runif(50,0,10);r3 <- runif(50,0,10)
apply(cbind(r1,r2,r3),1, FUN = function(x){ll1.wei(x,Ti[,1])})
table2 <- matrix(0,nrow = 8,ncol =7)
system.time({
  for(k in 1:8){
    optim.result=matrix(0,nrow = 1,ncol = 7)
    for (i in 1:50){
      fit=optim(c(r1[i],r2[i],r3[i]),ll1.wei,hessian=T,data=Ti[,k])
      fit1=optim(fit$par,ll1.wei,hessian=T,data=Ti[,k])
      while (abs(fit1$value-fit$value)>10^(-3)){
        fit <- fit1
        fit1<- optim(fit$par,ll1.wei,hessian=T,data=Ti[,k])
      }
      I=solve(fit1$hessian)
      sda <- sqrt(I[1,1]*10^(-4))
      sdb <- sqrt(I[2,2]*10^(-8))
      sds <- sqrt(I[3,3]*10^(4))
      t=cbind(fit1$par[1]*10^(-2),sda,
              fit1$par[2]*10^(-4),sdb,
              abs(fit1$par[3]*10^(2)),sds,
              (-fit1$value))
      optim.result=rbind(optim.result,t)
    }
    optim.result <- optim.result[-1,]
    print(optim.result)
    ind <- which.max(optim.result[,7])
    table2[k,] <- optim.result[ind,]
  }
})
colnames(table2) <- c("a","sda","b","sdb","sigma","sds","logL")
table2.wei <- table2[,c(1,3,5,7)]

#-------------------------------------------------------------------------------
#Figure 5
#-------------------------------------------------------------------------------

S.wei <- data.frame(S=table2[,c(1,3,5)],
                    group=c(1,1,1,3,3,3,5,5))
lma.wei <- lm(S.a~group,data=S.wei[,c(1,4)])
lmb.wei <- lm(S.b~group,data=S.wei[,c(2,4)])
lms.wei <- lm(S.sigma~group,data=S.wei[,c(3,4)])


par(mfrow=c(1,3))
plot(S.wei$group,S.wei$S.a,xlab = "Discharge Rate",ylab = "a",main="Weibull")
abline(lma.wei)

plot(S.wei$group,S.wei$S.b,xlab = "Discharge Rate",ylab = "b",main="Weibull")
abline(lmb.wei)

plot(S.wei$group,S.wei$S.sigma,xlab = "Discharge Rate",
     ylab = expression(beta),main="Weibull")
abline(lms.wei)
par(mfrow=c(1,1))
table3.wei <- rbind(summary(lma.wei)$coef[,1:2],summary(lmb.wei)$coef[,1:2],summary(lms.wei)$coef[,1:2])
row.names(table3.wei) <- c(expression(a_0),expression(a_1),
                           expression(b_0),expression(b_1),
                           expression(beta_0),expression(beta_1))

#------------------- use regression method to  sumulate table3 ------------------#

table3.wei<- cbind(table3.wei,rbind(confint(lma.wei),confint(lmb.wei),confint(lms.wei)))


#---------------------------------------- write output -------------------------------------------#
result <- list(log=table3.log,
               wei=table3.wei,
               nor=table3.nor)
library(data.table)
outputfile <- "result.csv" #output file name
sep <- "," #define the separator (related to format of the output file)
for(i in names(result)){
  fwrite(list(i), file=outputfile, sep=sep, append=T) #write names of the list elements
  ele <- result[[i]]
  fwrite(ele, file=outputfile, sep=sep, append=T, col.names=T) 
  fwrite(list(NA), file=outputfile, append=T) #add an empty row to separate elements
}

#----------------------------------------------------------------------------------------------------#

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
ind=which(result[,7]>624.5349)
result1 <- result[ind,]
ind1 <- which.max(result1[,7])
t3.wei <- matrix(c(result1[ind1,1]*10^(-2),result1[ind1,2]*10^(-3),result1[ind1,3]*10^(-4),
                       result1[ind1,4]*10^(-5),result1[ind1,5],result1[ind1,6],result1[ind1,7]),ncol=1)
tll.wei(0.091455547,
        0.001016102,
        0.000209264,
        1.15E-05,
        548.0750741,
        -49.92868114,Ti)

#table3
#---------------------
#          MLE
#---------------------
# a0   0.09150811
# a1   0.001019208
# b0   0.0002068806
# b1   1.163571e-05
# s0   3.155421e+02
# s1  -2.167785e+01
#logL  6.245349e+02
#---------------------
#---------------------
#      MLE(同論文 -> 四捨五入到小數點第二位)
#---------------------
# a0   9.15*10^(-2)
# a1   1.02*10^(-3)
# b0   2.07*10^(-4)
# b1   1.16*10^(-5)
# s0   3.15*10^2
# s1   -2.17*10
#logL   6.245349e+02
#---------------------




#---------------------------------------- write output -------------------------------------------#
result <- list(log=t3.log,
               wei=t3.wei,
               nor=t3.nor)
library(data.table)
outputfile <- "result.csv" #output file name
sep <- "," #define the separator (related to format of the output file)
for(i in names(result)){
  fwrite(list(i), file=outputfile, sep=sep, append=T) #write names of the list elements
  ele <- result[[i]]
  fwrite(ele, file=outputfile, sep=sep, append=T, col.names=T) 
  fwrite(list(NA), file=outputfile, append=T) #add an empty row to separate elements
}







### ----------------------------------------------------------------------------
### table 1
### ----------------------------------------------------------------------------
t1ll <- function(theta,data){ # ll is -loglikeihood, so here is also -logL
  a0 <- theta[1];a1 <- theta[2];b0 <- theta[3];b1 <- theta[4];s0 <- theta[5];s1 <- theta[6]
  aS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) a0+a1*x)
  bS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) b0+b1*x)
  sS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) s0+s1*x)
  t1 <- sum(apply(matrix(1:3,ncol = 1),1,function(x) ll(aS[1],bS[1],sS[1],data[,x])))
  t2 <- sum(apply(matrix(4:6,ncol = 1),1,function(x) ll(aS[2],bS[2],sS[2],data[,x])))
  t3 <- sum(apply(matrix(7:8,ncol = 1),1,function(x) ll(aS[3],bS[3],sS[3],data[,x])))
  return(c(-t1,-t2,-t3))
}

t1ll.log <- function(theta,data){ # ll is -loglikeihood, so here is also -logL
  a0 <- theta[1];a1 <- theta[2];b0 <- theta[3];b1 <- theta[4];s0 <- theta[5];s1 <- theta[6]
  aS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) a0+a1*x)
  bS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) b0+b1*x)
  sS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) s0+s1*x)
  t1 <- sum(apply(matrix(1:3,ncol = 1),1,function(x) ll.log(aS[1],bS[1],sS[1],data[,x])))
  t2 <- sum(apply(matrix(4:6,ncol = 1),1,function(x) ll.log(aS[2],bS[2],sS[2],data[,x])))
  t3 <- sum(apply(matrix(7:8,ncol = 1),1,function(x) ll.log(aS[3],bS[3],sS[3],data[,x])))
  return(c(-t1,-t2,-t3))
}

t1ll.wei<- function(theta,data){ # ll is -loglikeihood, so here is also -logL
  a0 <- theta[1];a1 <- theta[2];b0 <- theta[3];b1 <- theta[4];s0 <- theta[5];s1 <- theta[6]
  aS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) a0+a1*x)
  bS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) b0+b1*x)
  sS <- apply(matrix(c(1,3,5),ncol = 1),1,function(x) s0+s1*x)
  t1 <- sum(apply(matrix(1:3,ncol = 1),1,function(x) ll.wei(aS[1],bS[1],sS[1],data[,x])))
  t2 <- sum(apply(matrix(4:6,ncol = 1),1,function(x) ll.wei(aS[2],bS[2],sS[2],data[,x])))
  t3 <- sum(apply(matrix(7:8,ncol = 1),1,function(x) ll.wei(aS[3],bS[3],sS[3],data[,x])))
  return(c(-t1,-t2,-t3))
}

table1.LN <- rbind(t1ll.wei(as.vector(t3.wei),Ti), t1ll.log(as.vector(t3.log),Ti), t1ll(as.vector(t3.nor),Ti))
row.names(table1.LN) <- c("Weibull","Log-normal","Normal")
colnames(table1.LN) <- c("1C","3C","5C")

table1.LN 



## power-law 

ll <- function(a,b,s,data){  #-loglikelihood
  sumlt <- sum(log(data[-1]));sq <- s^2
  k <- c()
  for (s in 2:J){
    t<- (a/b*(exp(b*data[s])-exp(b*data[s-1]))-1)^2
    k <- append(k,t)
  }
  ans= -45/2*log(2*pi*sq)+45*log(a)+b*sumt-(sum(k)/(2*sq))
  -ans
}
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
apply(cbind(r1,r2,r3,r4,r5,r6),1, FUN = function(x){tll1(x,data=Ti)})
#apply(cbind(r1,r2,r3,r4,r5,r6),1, FUN = function(x){optim(x,tll1,hessian=F,data=Ti)})
start_time <- Sys.time()
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
end_time <- Sys.time()
end_time - start_time
optim.result <- optim.result[-1,]
ind <- which.max(optim.result[,7])
t3.nor <- optim.result[ind,]
t3.nor <- matrix(t3.nor,ncol = 1)
t3sim <- optim.result[ind,1:6]














