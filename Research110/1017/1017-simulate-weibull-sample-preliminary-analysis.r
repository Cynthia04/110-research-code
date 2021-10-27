t3 <- matrix(c(3,
               2,
               5,
               4,
               20,
               -5),ncol=2,byrow = T)
sk <- c(0.111,0.556,1)
est <- apply(t3,1,function(x) x[1]+sk*x[2])
colnames(est) <- c("a","b","s")
rownames(est) <- c("1C","3C","5C")
n1 <- 3 ;n2 <- 3 ;n3 <- 3
nk <- n1+n2+n3
est0 <- rbind(matrix(rep(est[1,],n1),ncol=3,byrow = T),
              matrix(rep(est[2,],n2),ncol=3,byrow = T),
              matrix(rep(est[3,],n3),ncol=3,byrow = T))
colnames(est0) <- c("a","b","s")
set.seed(1)
x <- apply(matrix(est0[,3],ncol=1),1,
           function(x){rweibull(45,shape=x,scale=(1/gamma(1/x+1)))})
x1 <- apply(x,2,function(s){cumsum(s)})
x2 <- matrix(0,ncol=3)
for(i in 1:nk){
  t <- cbind (rep(est0[i,1],45),rep(est0[i,2],45),x1[,i])
  x2 <- rbind(x2,t)
}
x2 <- x2[-1,]
Ti<- apply(x2,1,function(x){log(x[2]*x[3]/x[1]+1)/x[2]})
Ti <- matrix(Ti,ncol=nk)
Ti <- rbind(rep(0,nk),Ti)
rm(i,t)


#----------------------------------------------

# Preliminary analysis

#----------------------------------------------
# 1. plot 

par(mfrow=c(1,2))
matplot(Ti,type="b",pch=1:3,col=1:3,
        xlab = "45 event",ylab = "time",
        main="arrival time of the 3 sample under normal use in the original process",
        cex.lab=2, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
asd <- Ti[2:46,]-Ti[1:45,]
matplot(asd,
        xlab = "45 event",ylab = "time",
        main="Interarrival times of the 3 sample under normal use in the original process",
        pch=1,col="black",
        cex.lab=2, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

par(mfrow=c(1,1))

#----------------------------------------------
f <- function(a,b,data){
  (a/b)*(exp(b*data[s])-exp(b*data[s-1]))}



sk <- c(0.111,0.556,1)
estresult <- apply(matrix(t3.wei[1:6,1],ncol=2,byrow = T),1,function(x) x[1]+sk*x[2])
est0result <- rbind(matrix(rep(estresult[1,],n1),ncol=3,byrow = T),
                    matrix(rep(estresult[2,],n2),ncol=3,byrow = T),
                    matrix(rep(estresult[3,],n3),ncol=3,byrow = T))

Xi <- apply(cbind(est0result[,1:2],1:9),1,function(x){x[1]/x[2]*(exp(x[2]*Ti[2:46,x[3]])-exp(x[2]*Ti[1:45,x[3]]))})
matplot(Xi,
        xlab = "45 event",ylab = "time",
        main="Interarrival times of the 9 sample in the transformed process",
        pch=c(rep(1,n1),rep(2,n2),rep(4,n3)),col="black")

est0result[,3]


#----------------------------------------------

#2. estimate the parameter of each sample to see the relationship between each sample

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
est0
ll.wei(est[1,1],est[1,2],est[1,3],Ti[,1])
# ll0 <- function(theta,data){
#   ll(theta[1],theta[2],theta[3],data)
# }

# ll1.wei <- function(theta,data){
#   ll(theta[1]*10^(-2),theta[2]*10^(-4),theta[3]*10^(-3),data)
# }
# 確認variance
# fit0 <- optim(c(9.219*10^(-2),2.182*10^(-4),2.339*10^(-3)),ll0,hessian=T,data=Ti[,1]);fit0
# fit1 <- optim(c(9.219,2.182,2.339),ll1.wei,hessian=T,data=Ti[,1]);fit1
# 
# 1/fit0$hessian[1,1]
# 10^(-4)*1/fit1$hessian[1,1]
# solve(fit0$hessian)[1,1]
# 10^(-4)*solve(fit1$hessian)[1,1]

ll1.wei <- function(theta,data){
  ll.wei(theta[1]*10^(-1),theta[2]*10^(-1),theta[3]*10^(-1),data)
}
ll1.wei(c(est[1,1]*10,est[1,2]*10,est[1,3]*10),Ti[,1])

set.seed(109225017);r1 <- runif(50,0,100);r2 <- runif(50,0,100);r3 <- runif(50,0,100)
apply(cbind(r1,r2,r3),1, FUN = function(x){ll1.wei(x,Ti[,1])})
#apply(cbind(r1,r2,r3),1, FUN = function(x){optim(x,ll1.wei,hessian=T,data=Ti[,1])})#確定無Inf or NA

table2 <- matrix(0,nrow = 9,ncol =7)
system.time({
  for(k in 1:9){
    optim.result=matrix(0,nrow = 1,ncol = 7)
    for (i in 1:50){
      fit=optim(c(r1[i],r2[i],r3[i]),ll1.wei,hessian=T,data=Ti[,k])
      fit1=optim(fit$par,ll1.wei,hessian=T,data=Ti[,k])
      while (abs(fit1$value-fit$value)>10^(-3)){
        fit <- fit1
        fit1<- optim(fit$par,ll1.wei,hessian=T,data=Ti[,k])
      }
      I=solve(fit1$hessian)
      sda <- sqrt(I[1,1]*10^(-2))
      sdb <- sqrt(I[2,2]*10^(-2))
      sds <- sqrt(I[3,3]*10^(-2))
      t=cbind(fit1$par[1]*10^(-1),sda,
              fit1$par[2]*10^(-1),sdb,
              fit1$par[3]*10^(-1),sds,
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
plotCI(1:9, est0[,1], li=ci[,1], ui=ci[,2],
       xlab="sample",ylab=expression(hat(a)),
       main="95% Confidence Interval of a")
plotCI(1:9, est0[,2], li=ci[,3], ui=ci[,4],
       xlab="sample",ylab=expression(hat(b)),
       main="95% Confidence Interval of b")
plotCI(1:9, est0[,3], li=ci[,5], ui=ci[,6],
       xlab="sample",ylab=expression(hat(beta)),
       main=expression(paste("95% Confidence Interval of ",beta)))
par(mfrow=c(1,1))
#-------------- calculate coverage probability & average length -----------------#

# Above is one sample of simulation
# We need to simulate large sample
# ---> 1017 simulate weibull sample cpal.r 

#-----------------------------------------------------------------------------------
#Figure 5

S.abs <- data.frame(S=table2[,c(1,3,5)],
                    group=c(rep(sk[1],3),rep(sk[2],3),rep(sk[3],3)))
lma <- lm(S.a~group,data=S.abs[,c(1,4)])
lmb <- lm(S.b~group,data=S.abs[,c(2,4)])
lms <- lm(S.sigma~group,data=S.abs[,c(3,4)])


par(mfrow=c(1,3))
plot(S.abs$group,S.abs$S.a,xlab = "Discharge Rate (standardization)",
     ylab = expression(hat(a)),main="Weibull",
     cex.lab=1.5, cex.axis=1, cex.main=2, cex.sub=1)
abline(lma)

plot(S.abs$group,S.abs$S.b,xlab = "Discharge Rate (standardization)",
     ylab = expression(hat(b)),main="Weibull",
     cex.lab=1.5, cex.axis=1, cex.main=2, cex.sub=1)
abline(lmb)

plot(S.abs$group,S.abs$S.sigma,xlab = "Discharge Rate (standardization)",
     ylab = expression(hat(beta)),main="Weibull",
     cex.lab=1.5, cex.axis=1, cex.main=2, cex.sub=1)
abline(lms)
par(mfrow=c(1,1))

summary(lma)
summary(lmb)
summary(lms)
table3 <- rbind(summary(lma)$coef[,1:2],summary(lmb)$coef[,1:2],summary(lms)$coef[,1:2])
row.names(table3) <- c(expression(a_0),expression(a_1),
                       expression(b_0),expression(b_1),
                       expression(beta_0),expression(beta_1))
table3



result <- list(table2=table2,
               ci=ci,
               table3=table3)
library(data.table)
outputfile <- "result.csv" #output file name
sep <- "," #define the separator (related to format of the output file)
for(i in names(result)){
  fwrite(list(i), file=outputfile, sep=sep, append=T) #write names of the list elements
  ele <- result[[i]]
  fwrite(ele, file=outputfile, sep=sep, append=T, col.names=T) 
  fwrite(list(NA), file=outputfile, append=T) #add an empty row to separate elements
}
