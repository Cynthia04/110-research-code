#---------------------
#      MLE(同論文 -> 四捨五入到小數點第二位)
#---------------------
# a0   9.15*10^(-2)
# a1   1.01*10^(-3)
# b0   2.07*10^(-4)
# b1   1.19*10^(-5)
# s0   2.76*10^(-3)
# s1   4.24*10^(-4)
#---------------------
#將以上結果視為真實參數 進行有母數拔靴法 生出 500 bootstrap samples
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
#---------------------------------------------------------------#
# 4.4 Confidence intervals for theta and EOP
#---------------------------------------------------------------#


f <- function(g){
  #likelihood fn 
  ll <- function(a,b,s,data){  #-loglikelihood
    sumt <- sum(data);J=length(data);sq <- s^2
    k <- c()
    for (s in 2:J){
      t<- (a/b*(exp(b*data[s])-exp(b*data[s-1]))-1)^2
      k <- append(k,t)
    }
    ans= -g/2*log(2*pi*sq)+g*log(a)+b*sumt-(sum(k)/(2*sq))
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
  x <- apply(matrix(est[,3],ncol = 1),1,
             function(x){c(rnorm(g,mean=1,sd=x),
                           rnorm(g,mean=1,sd=x),
                           rnorm(g,mean=1,sd=x))})
  x <- cbind(matrix(x[,1],ncol = 3),matrix(x[,2],ncol = 3),matrix(x[,3],ncol = 3))
  x <- x[,-9]
  x1 <- apply(x,2,function(s){cumsum(s)})
  x2 <- matrix(0,ncol=3)
  for(i in 1:8){
    t <- cbind (rep(est0[i,1],g),rep(est0[i,2],g),x1[,i])
    x2 <- rbind(x2,t)
  }
  x2 <- x2[-1,]
  Ti<- apply(x2,1,function(x){log(x[2]*x[3]/x[1]+1)/x[2]})
  Ti <- matrix(Ti,ncol=8)
  Ti <- rbind(rep(0,8),Ti)
  #estimate
  r1 <- runif(10,0,100);r2 <- runif(10,0,100);r3 <- runif(10,0,100);r4 <- runif(10,0,100);r5 <- runif(10,0,100);r6 <- runif(10,0,100)
  optim.result=matrix(0,nrow = 1,ncol = 7)
  for (i in 1:10){
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
  t3sim <- optim.result[ind,1:6]
  theta <- as.numeric(t3sim)
  EOP <- matrix(0,ncol =1,nrow = 1)
  S0 <- c()
  for (i in c(1,3,5)){
    t <- t3sim[i]+0.5*t3sim[i+1]
    S0 <- append(S0,t)
  }
  EZ <- function(i){
    a <- S0[1];b <- S0[2];s <- S0[3]
    k <- a/b
    1/b*(log((i+k)/(i-1+k))+s^2/2*((i-1)/(i-1+k)^2-i/(i+k)^2))
  }
  EOP=1; ans <- EZ(EOP)
  while( ans>8){
    EOP = EOP+1 
    ans = EZ(EOP)
  }
  return(c(theta,EOP))
}

start_time <- Sys.time()
result <- replicate(500,f(45))
end_time <- Sys.time()
end_time - start_time









