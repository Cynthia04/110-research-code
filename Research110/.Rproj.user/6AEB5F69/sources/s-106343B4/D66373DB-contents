#-------------------------------------------------------------------------------------------------
# 設定參數 ---------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------------------------
# 生資料 -----------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
set.seed(1) #fix one data 
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

#-------------------------------------------------------------------------------------------------
# mle 估計 ---------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
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


tll.wei<- function(a0,a1,b0,b1,s0,s1,data){ # ll is -loglikeihood, so here is also -logL
  aS <- apply(matrix(sk,ncol = 1),1,function(x) a0+a1*x)
  bS <- apply(matrix(sk,ncol = 1),1,function(x) b0+b1*x)
  sS <- apply(matrix(sk,ncol = 1),1,function(x) s0+s1*x)
  t1 <- sum(apply(matrix(1:n1,ncol = 1),1,function(x) ll.wei(aS[1],bS[1],sS[1],data[,x])))
  t2 <- sum(apply(matrix((n1+1):(n1+n2),ncol = 1),1,function(x) ll.wei(aS[2],bS[2],sS[2],data[,x])))
  t3 <- sum(apply(matrix((n1+n2+1):nk,ncol = 1),1,function(x) ll.wei(aS[3],bS[3],sS[3],data[,x])))
  ans <- t1+t2+t3
  ans
}
tll1.wei <- function(theta,data){
  tll.wei(theta[1]*10^(-1),theta[2]*10^(-1),theta[3]*10^(-1),theta[4]*10^(-1),theta[5]*10^(-1),theta[6]*10^(-1),data)
}
Opt<-function(x){
  iter=x[7]
  initial=c(x[1],x[2],x[3],x[4],x[5],x[6])
  error=1
  before=0
  while(abs(error)-1e-15>0){
    step1=try(optim(initial,tll1.wei,control=list(factr=1e-15,maxit=10^9),hessian=T,data=Ti),silent = T)
    if(step1$convergence==0){
      initial <- step1$par
      error <- step1$value-before
      before <- step1$value
    }
  }
  parameter <- c(step1$par[1]*10^(-1),step1$par[2]*10^(-1),step1$par[3]*10^(-1),step1$par[4]*10^(-1),step1$par[5]*10^(-1),step1$par[6]*10^(-1))
  I=solve(step1$hessian)
  sd <- c(sqrt(I[1,1]*10^(-2)),sqrt(I[2,2]*10^(-2)),sqrt(I[3,3]*10^(-2)),
          sqrt(I[4,4]*10^(-2)),sqrt(I[5,5]*10^(-2)),sqrt(I[6,6]*10^(-2)))
  return(c(parameter,-step1$value,sd)) #logL
}
N=50
rn <- matrix(0,ncol=7,nrow=N)
set.seed(1);rn[,1]<-runif(N,0,100);rn[,2]<-runif(N,0,100);rn[,3]<-runif(N,0,100)
set.seed(1);rn[,4]<-runif(N,0,100);rn[,5]<-runif(N,0,100);rn[,6]<-runif(N,0,100)
rn[,7]<-1:N
#Opt(as.matrix(rn)[1,])
library(snow)
library(pbapply)
library(parallel)
cl = makeCluster(4,type='SOCK') 
clusterExport(cl, c("n1","n2","n3","nk","sk","rn","Opt","ll.wei","tll.wei","tll1.wei","Ti")) 
result=pbapply(as.matrix(rn),1,function(x){Opt(x)},cl = cl)
stopCluster(cl)
result <- t(result)
ind1 <- which.max(result[,7])
t3.wei <- matrix(result[ind1,],ncol=1)

#-------------------------------------------------------------------------------------------------
# H-M algorithm ---------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
tll.wei<- function(theta,data){ # ll is -loglikeihood, -ans -> here is logL
  a0 <- theta[1];a1 <- theta[2];b0 <- theta[3];b1 <- theta[4];s0 <- theta[5];s1 <- theta[6]
  aS <- apply(matrix(sk,ncol = 1),1,function(x) a0+a1*x)
  bS <- apply(matrix(sk,ncol = 1),1,function(x) b0+b1*x)
  sS <- apply(matrix(sk,ncol = 1),1,function(x) s0+s1*x)
  t1 <- sum(apply(matrix(1:n1,ncol = 1),1,function(x) ll.wei(aS[1],bS[1],sS[1],data[,x])))
  t2 <- sum(apply(matrix((n1+1):(n1+n2),ncol = 1),1,function(x) ll.wei(aS[2],bS[2],sS[2],data[,x])))
  t3 <- sum(apply(matrix((n1+n2+1):nk,ncol = 1),1,function(x) ll.wei(aS[3],bS[3],sS[3],data[,x])))
  ans <- t1+t2+t3
  -ans
}

m=5*10^4;N=5000;k=10
a00 <- 1;a10 <- 1;b00 <- 1;b10 <- 1;s00 <-13;s10 <- -2
sigma_proposal=c(0.05,0.05,0.05,0.05,2,1);sigma_prior=c(0.1,0.5,0.1,0.1,3,1) 
a0 <- c();a1 <- c();b0 <- c();b1 <- c();s0 <- c();s1 <- c()
alpha1 <- c(); alpha2 <- c(); alpha3 <- c(); alpha4 <- c(); alpha5 <- c(); alpha6 <- c()
u <- runif(m+N*k)
start <- Sys.time()
for(i in 1:(m+N*k)){
  #a0
  repeat{
    z1 <- rnorm(1,a00,sigma_proposal[1])
    if(z1>0){break}
  }
  p1 <- exp(-(z1-t3.wei[1,1])^2/(2*sigma_prior[1]^2))
  p0 <- exp(-(a00-t3.wei[1,1])^2/(2*sigma_prior[1]^2))
  alpha1[i] <- alpha <- min(exp(tll.wei(c(z1,a10,b00,b10,s00,s10),Ti)-tll.wei(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha>u[i]){
    a0[i] <- z1
  }else{
    a0[i] <- a00
  }
  a00 <- a0[i]
  #a1
  repeat{
    z2 <- rnorm(1,a10,sigma_proposal[2])
    if(z2>0){break}
  }
  p1 <- exp(-(z2-t3.wei[2,1])^2/(2*sigma_prior[2]^2))
  p0 <- exp(-(a10-t3.wei[2,1])^2/(2*sigma_prior[2]^2))
  alpha2[i] <- alpha <- min(exp(tll.wei(c(a00,z2,b00,b10,s00,s10),Ti)-tll.wei(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha>u[i]){
    a1[i] <- z2
  }else{
    a1[i] <- a10
  }
  a10 <- a1[i]
  #b0
  repeat{
    z3 <- rnorm(1,b00,sigma_proposal[3])
    if(z3>0){break}
  }
  p1 <- -(z3-t3.wei[3,1])^2/(2*sigma_prior[3]^2)
  p0 <- -(b00-t3.wei[3,1])^2/(2*sigma_prior[3]^2)
  alpha3[i] <- alpha <- min(exp(tll.wei(c(a00,a10,z3,b10,s00,s10),Ti)-tll.wei(c(a00,a10,b00,b10,s00,s10),Ti))*exp(p1-p0),1)
  if(alpha>u[i]){
    b0[i] <- z3
  }else{
    b0[i] <- b00
  }
  b00 <- b0[i]
  #b1
  repeat{
    z4 <- rnorm(1,b10,sigma_proposal[4])
    if(z4>0){break}
  }
  p1 <- exp(-(z4-t3.wei[4,1])^2/(2*sigma_prior[4]^2))
  p0 <- exp(-(b10-t3.wei[4,1])^2/(2*sigma_prior[4]^2))
  alpha4[i] <- alpha <- min(exp(tll.wei(c(a00,a10,b00,z4,s00,s10),Ti)-tll.wei(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha>u[i]){
    b1[i] <- z4
  }else{
    b1[i] <- b10
  }
  b10 <- b1[i]
  #beta0
  repeat{
    z5 <- rnorm(1,s00,sigma_proposal[5])
    if(z5>0){break}
  }
  p1 <- -(z5-t3.wei[5,1])^2/(2*sigma_prior[5]^2)
  p0 <- -(s00-t3.wei[5,1])^2/(2*sigma_prior[5]^2)
  alpha5[i] <- alpha <- min(exp(tll.wei(c(a00,a10,b00,b10,z5,s10),Ti)-tll.wei(c(a00,a10,b00,b10,s00,s10),Ti))*exp(p1-p0),1)
  if(alpha>u[i]){
    s0[i] <- z5
  }else{
    s0[i] <- s00
  }
  s00 <- s0[i]
  #beta1
  repeat{
    z6 <- rnorm(1,s10,sigma_proposal[6])
    if(z6>-7){break}
  }
  p1 <- -(z6-t3.wei[6,1])^2/(2*sigma_prior[6]^2)
  p0 <- -(s10-t3.wei[6,1])^2/(2*sigma_prior[6]^2)
  alpha6[i] <- alpha <- min(exp(tll.wei(c(a00,a10,b00,b10,s00,z6),Ti)-tll.wei(c(a00,a10,b00,b10,s00,s10),Ti))*exp(p1-p0),1)
  if(alpha>u[i]){
    s1[i] <- z6
  }else{
    s1[i] <- s10
  }
  s10 <- s1[i]
}
end <- Sys.time()
end-start

#----------------------------------------------------------------------------------------------
# H-M algorithm 
# iteraction times 1e+05 | burn in size 5e+04 | Thin 10

# normal proposal density
# initial value true value (3,2,5,4,20,-5) | sd of  proposal (3,3,3,3,3,3)

# normal prior 
#mean of prior mle = (t3.wei[1:6,1])
# sd of prior (0.1,0.5,0.1,0.1,2,3) 
# result (running time  32.87364 mins)
#----------------------------------------------------------------------------------------------
mcmc <- function(x){return(c(mean(x),sd(x),quantile(x, probs = c(0.025, 0.975))))}
mcmc_table <- t(apply(rbind(a0[seq(m+k,m+N*k,k)],a1[seq(m+k,m+N*k,k)],b0[seq(m+k,m+N*k,k)],
                            b1[seq(m+k,m+N*k,k)],s0[seq(m+k,m+N*k,k)],s1[seq(m+k,m+N*k,k)]),1,mcmc))

mcmc_table <- cbind(initial=c(3,2,5,4,20,-5),
                    sigma_proposal,
                    sigma_prior,
                    Acceptance=apply(rbind(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6),1,mean),
                    True=c(3,2,5,4,20,-5),mean=mcmc_table[,1],mle=t3.wei[1:6,1],
                    sd=mcmc_table[,2],"sd of mle"=t3.wei[8:13,1],mcmc_table[,-c(1:2)])
mcmc_table
write.csv(mcmc_table,"mcmc.csv")

par(mfrow=c(3,1))
plot(a0,type = "l",xlab = "iteration times (1e+05)",ylab="a0",main="Trace of MCMC sample of a0",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(a1,type = "l",xlab = "iteration times (1e+05)",ylab="a1",main="Trace of MCMC sample of a1",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(b0,type = "l",xlab = "iteration times (1e+05)",ylab="b0",main="Trace of MCMC sample of b0",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(b1,type = "l",xlab = "iteration times (1e+05)",ylab="b1",main="Trace of MCMC sample of b1",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(s0,type = "l",xlab = "iteration times (1e+05)",ylab="beta0",main="Trace of MCMC sample of beta0",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(s1,type = "l",xlab = "iteration times (1e+05)",ylab="beta1",main="Trace of MCMC sample beta1",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
par(mfrow=c(1,1))

par(mfrow=c(3,1))
plot(a0[seq(m+k,m+N*k,k)],type = "l",xlab = "iteration times (5000)",
     ylab="a0",main="Trace of posterior sample of a0",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(a1[seq(m+k,m+N*k,k)],type = "l",xlab = "iteration times (5000)",
     ylab="a1",main="Trace of posterior sample of a1",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(b0[seq(m+k,m+N*k,k)],type = "l",xlab = "iteration times (5000)",
     ylab="b0",main="Trace of posterior sample of b0",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(b1[seq(m+k,m+N*k,k)],type = "l",xlab = "iteration times (5000)",
     ylab="b1",main="Trace of posterior sample of b1",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(s0[seq(m+k,m+N*k,k)],type = "l",xlab = "iteration times (5000)",
     ylab="beta0",main="Trace of posterior sample of beta0",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(s1[seq(m+k,m+N*k,k)],type = "l",xlab = "iteration times (5000)",
     ylab="beta1",main="Trace of posterior sample of beta1",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
par(mfrow=c(1,1))

par(mfrow=c(3,1))
plot(density(a0[seq(m+k,m+N*k,k)]),main="posterior density of a0",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(density(a1[seq(m+k,m+N*k,k)]),main="posterior density of a1",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(density(b0[seq(m+k,m+N*k,k)]),main="posterior density of b0",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(density(b1[seq(m+k,m+N*k,k)]),main="posterior density of b1",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(density(s0[seq(m+k,m+N*k,k)]),main="posterior density of beta0",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(density(s1[seq(m+k,m+N*k,k)]),main="posterior density of beta1",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
par(mfrow=c(1,1))
