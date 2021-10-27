#未加速 (正常應力下)去做
#即僅估計截距項參數: a0,b0,beta0
library(snow)
library(pbapply)
library(parallel)
t3 <- matrix(c(3,
               2,
               5,
               4,
               20,
               -5),ncol=2,byrow = T)
sk <- 0
est <- matrix(apply(t3,1,function(x) x[1]+sk*x[2]),nrow=1)
colnames(est) <- c("a","b","s")
rownames(est) <- c("0.5C")
n0=3
est0 <- matrix(rep(est[1,],n0),ncol=3,byrow = T)
colnames(est0) <- c("a","b","s")


#------------------------------------------------------------------------
# estimate parameter of general likelihood
# likelihood function of one sample
#------------------------------------------------------------------------
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


tll.wei<- function(a,b,s,data){ # ll is -loglikeihood, so here is also -logL
  t1 <- sum(apply(matrix(1:n0,ncol = 1),1,function(x) ll.wei(a,b,s,data[,x])))
  ans <- t1
  ans
}
tll1.wei <- function(theta,data){
  tll.wei(theta[1]*10^(-1),theta[2]*10^(-1),theta[3]*10^(-1),data)
}
Opt<-function(x){
  iter=x[4]
  initial=c(x[1],x[2],x[3])
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
  parameter <- c(step1$par[1]*10^(-1),step1$par[2]*10^(-1),step1$par[3]*10^(-1))
  I=solve(step1$hessian)
  sd <- c(sqrt(I[1,1]*10^(-2)),sqrt(I[2,2]*10^(-2)),sqrt(I[3,3]*10^(-2)))
  return(c(parameter,-step1$value,sd)) #logL
}
start <- Sys.time()
Opt(as.matrix(rn)[1,])
end<- Sys.time()
end-start

#---  method = "BFGS"
#---  跑一次要18分鐘相比於default 只要3秒



################################################################################################################################
par0 <- matrix(0,ncol=1,nrow = 7)
for( seed in 128:250){
  set.seed(seed)
  x <- apply(matrix(est0[,3],ncol=1),1,
             function(x){rweibull(45,shape=x,scale=(1/gamma(1/x+1)))})
  x1 <- apply(x,2,function(s){cumsum(s)})
  x2 <- matrix(0,ncol=3)
  for(i in 1:n0){
    t <- cbind (rep(est0[i,1],45),rep(est0[i,2],45),x1[,i])
    x2 <- rbind(x2,t)
  }
  x2 <- x2[-1,]
  Ti<- apply(x2,1,function(x){log(x[2]*x[3]/x[1]+1)/x[2]})
  Ti <- matrix(Ti,ncol=n0)
  Ti <- rbind(rep(0,n0),Ti)
  rm(i,t)
  N=50
  rn <- matrix(0,ncol=4,nrow=N)
  set.seed(seed);rn[,1]<-runif(N,0,100);rn[,2]<-runif(N,0,100);rn[,3]<-runif(N,0,100)
  rn[,4]<-1:N
  #Opt(as.matrix(rn)[1,])
  result=apply(as.matrix(rn),1,function(x){Opt(x)})
  result <- t(result)
  ind1 <- which.max(result[,4])
  t3.wei <- matrix(result[ind1,],ncol=1)
  par0 <- cbind(par0,t3.wei)
  print(seed)
}

finalpar0 <- cbind(finalpar0,par0[,-1])


# seed: 1:16、18:102、104:106、108-114、116:124、126

16+85









