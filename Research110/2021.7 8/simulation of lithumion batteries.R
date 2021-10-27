# data <-  data.frame(b1 = c(10.858,10.820,10.773,10.767,9.867,9.842),
#                     b2 = c(10.805,10.777,10.722,10.715,9.815,9.755),
#                     b3 = c(10.840,10.791,10.758,10.745,9.800,9.785),
#                     b4 = c(10.608,10.555,10.517,10.515,9.568,9.544),
#                     b5 = c(10.524,10.510,10.453,10.453,9.516,9.444),
#                     b6 = c(10.582,10.530,10.505,10.502,9.531,9.518),
#                     b7 = c(10.318,10.247,10.221,10.219,9.141,9.179),
#                     b8 = c(10.345,10.264,10.230,10.195,9.290,9.275))

#--------------- Use the normal estimate of table 2 to simulate data -----------#
est <- data.frame(a=c(9.219,9.269,9.235,9.422,9.503,9.434,9.649,9.638),
                  b=c(2.182,2.215,2.258,2.435,2.394,2.449,2.771,2.583),
                  sigma=c(2.339,2.080,1.782,3.357,2.518,2.541,3.944,4.053))
#----------------- Simulate data m_jk=45 (n_1,n_2,n_3)=(3,3,2) ----------------#
set.seed(1109225017)
x <- apply(matrix(c(est[,3]*10^(-3)),ncol = 1),1,function(x){rnorm(45,mean=1,sd=x)})
x <- apply(x,2,function(s){cumsum(s)})
x1 <- matrix(0,ncol=3)
for(i in 1:8){
  t <- cbind (rep(est[i,1],45),rep(est[i,2],45),x[,i])
  x1 <- rbind(x1,t)
}
x1 <- x1[-1,]
Ti<- apply(x1,1,function(x){log(x[2]*x[3]*10^(-2)/x[1]+1)/(x[2]*10^(-4))})
Ti <- matrix(Ti,ncol=8)
Ti <- rbind(rep(0,8),Ti)
rm(i,t)

#----------------------------- plot Figure 2 ------------------------------#

simulationdata1 <- Ti[2:46,]-Ti[1:45,]
asd <- simulationdata1
for(i in 1:8){
  asd[1:4,i] <- data[1:4,i]
}

for(i in 1:8){
  asd[44:45,i] <- data[5:6,i]
}

simud1.1 <- data.frame(y=as.vector(asd[,1:3]),
                       x=rep(1:45,3))
simud1.3 <- data.frame(y=as.vector(asd[,4:6]),
                       x=rep(1:45,3))
simud1.5 <- data.frame(y=as.vector(asd[,7:8]),
                       x=rep(1:45,2))
lm1.1 <- lm(y~x,data=simud1.1)
lm1.3 <- lm(y~x,data=simud1.3)
lm1.5 <- lm(y~x,data=simud1.5)

plot(1:45,asd[,1],xlim = c(0,45),ylim = c(9.0,11.0),xlab = "g-cycle",ylab = expression(CR[g]))
points(1:45,asd[,2])
points(1:45,asd[,3])
points(1:45,asd[,4], pch =2)
points(1:45,asd[,5], pch =2)
points(1:45,asd[,6], pch =2)
points(1:45,asd[,7], pch =4)
points(1:45,asd[,8], pch =4)

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
ll(9.219*10^(-2),2.182*10^(-4),2.339*10^(-3),Ti[,1])
# ll0 <- function(theta,data){
#   ll(theta[1],theta[2],theta[3],data)
# }

# ll1 <- function(theta,data){
#   ll(theta[1]*10^(-2),theta[2]*10^(-4),theta[3]*10^(-3),data)
# }
# 確認variance
# fit0 <- optim(c(9.219*10^(-2),2.182*10^(-4),2.339*10^(-3)),ll0,hessian=T,data=Ti[,1]);fit0
# fit1 <- optim(c(9.219,2.182,2.339),ll1,hessian=T,data=Ti[,1]);fit1
# 
# 1/fit0$hessian[1,1]
# 10^(-4)*1/fit1$hessian[1,1]
# solve(fit0$hessian)[1,1]
# 10^(-4)*solve(fit1$hessian)[1,1]

ll1 <- function(theta,data){
  ll(theta[1]*10^(-3),theta[2]*10^(-3),theta[3]*10^(-3),data)
}
ll1(c(9.219*10,2.182*10^(-1),2.339),Ti[,1])

set.seed(109225017);r1 <- runif(50,0,100);r2 <- runif(50,0,100);r3 <- runif(50,0,100)
apply(cbind(r1,r2,r3),1, FUN = function(x){ll1(x,Ti[,1])})
#apply(cbind(r1,r2,r3),1, FUN = function(x){optim(x,ll1,hessian=T,data=Ti[,1])})#確定無Inf or NA

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
table2

#ci
ci <- cbind(t(apply(table2[,1:2],1,function(x){x[1]+c(-1,1)*1.96*x[2]})),
            t(apply(table2[,3:4],1,function(x){x[1]+c(-1,1)*1.96*x[2]})),
            t(apply(table2[,5:6],1,function(x){x[1]+c(-1,1)*1.96*x[2]})))
as.numeric(ci)
colnames(ci) <- c("al","au","bl","bu","sl","su")

#plot of CI
require(plotrix)
par(mfrow=c(1,3))
plotCI(1:8, est[,1]*10^(-2), li=ci[,1], ui=ci[,2],
       xlab="sample",ylab=expression(hat(a)),
       main="95% Confidence Interval of a")
plotCI(1:8, est[,2]*10^(-4), li=ci[,3], ui=ci[,4],
       xlab="sample",ylab=expression(hat(b)),
       main="95% Confidence Interval of b")
plotCI(1:8, est[,3]*10^(-3), li=ci[,5], ui=ci[,6],
       xlab="sample",ylab=expression(hat(sigma)),
       main=expression(paste("95% Confidence Interval of ",sigma)))
par(mfrow=c(1,1))

#-------------- calculate coverage probability & average length -----------------#

# --- >>> see table2 cpal.r

#-----------------------------------------------------------------------------------

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
     ylab = expression(hat(sigma)),main="normal")
abline(lms)
par(mfrow=c(1,1))


table3 <- rbind(summary(lma)$coef[,1:2],summary(lmb)$coef[,1:2],summary(lms)$coef[,1:2])
row.names(table3) <- c(expression(a_0),expression(a_1),
                       expression(b_0),expression(b_1),
                       expression(sigma_0),expression(sigma_1))

#------------------- use regression method to  sumulate table3 ------------------#

table3 <- cbind(table3,rbind(confint(lma),confint(lmb),confint(lms)))
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
optim.result[ind,]
#---------------------------------------- write output -------------------------------------------#
result <- list(table2=table2,
               ci=ci,
               table3=table3,
               optim.result=optim.result)
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

t3sim <- optim.result[ind,1:6]

t3sim[1];t3sim[2];t3sim[3];t3sim[4];t3sim[5];t3sim[6]
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














