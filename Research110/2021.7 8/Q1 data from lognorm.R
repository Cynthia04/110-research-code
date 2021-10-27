#------
#0831
#Q1
#------
#lognormal
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
           function(x){rlnorm(45,meanlog=-x^2/2,sdlog=x^2)})
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
rm(i,t)

matplot(Ti,type="b",pch=1:8,col=1:8,
        xlab = "g-cycle",ylab = expression(CR[g]),
        main="Recurrent event proceses for the battery simulation data")
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

ll.log(9.219*10^(-2),2.182*10^(-4),2.34*10^(-3),Ti[,1])


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
colnames(table2) <- c("a","sda","b","sdb","sigma","sds","logL")
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
ll.wei(9.232*10^(-2),2.144*10^(-4),3.964*10^(2),Ti[,1])
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
              abs(fit1$par[3]*10^(2)),
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
table2.wei <- table2
sum(table2.wei[,4])

sum(table2.nor[,4])

sum(table2.log[,4])

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



#-----------------------------------------------------------------------------
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
    step1=try(optim(initial,tll1,control=list(factr=1e-15,maxit=10^9),data=Ti),silent = T)
    if(step1$convergence==0){
      initial <- step1$par
      error <- step1$value-before
      before <- step1$value
    }
  }
  return(c(step1$par,-step1$value)) 
}

library(snow)
library(pbapply)
library(parallel)
cl = makeCluster(4,type='SOCK') 
clusterExport(cl, c("rn","Opt","ll","tll","tll1","Ti")) 
result=pbapply(as.matrix(rn),1,function(x){Opt(x)},cl = cl)
stopCluster(cl)







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

library(snow)
library(pbapply)
library(parallel)
cpu.cores <- detectCores() 
cl = makeCluster(4,type='SOCK') 
clusterExport(cl, c("rn","Opt","ll.wei","tll.wei","tll1.wei","Ti")) 
result=pbapply(as.matrix(rn),1,function(x){Opt(x)},cl = cl)
stopCluster(cl)



result <- t(result)
write.csv(result,"log_data_wei_model.csv")
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


tll1
tll(theta[1],theta[2],theta[3],theta[4],theta[5],theta[6],Ti)


tll.log(theta[1],theta[2],theta[3],theta[4],theta[5],theta[6],Ti)



theta <- as.vector(t(t3))
