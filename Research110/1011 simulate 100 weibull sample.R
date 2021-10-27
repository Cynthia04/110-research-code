library(snow)
library(pbapply)
library(parallel)
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
start <- Sys.time()
Opt(as.matrix(rn)[1,])
end<- Sys.time()
end-start

#---  method = "BFGS"
#---  跑一次要18分鐘相比於default 只要3秒



################################################################################################################################

par <- matrix(0,ncol=1,nrow = 13)
for( seed in 1:100){
  seed=1
  set.seed(seed)
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
  
  N=50
  rn <- matrix(0,ncol=7,nrow=N)
  set.seed(seed);rn[,1]<-runif(N,0,100);rn[,2]<-runif(N,0,100);rn[,3]<-runif(N,0,100)
  set.seed(seed);rn[,4]<-runif(N,0,100);rn[,5]<-runif(N,0,100);rn[,6]<-runif(N,0,100)
  rn[,7]<-1:N
  #Opt(as.matrix(rn)[1,])
  cl = makeCluster(4,type='SOCK') 
  clusterExport(cl, c("n1","n2","n3","nk","sk","rn","Opt","ll.wei","tll.wei","tll1.wei","Ti")) 
  result=pbapply(as.matrix(rn),1,function(x){Opt(x)},cl = cl)
  stopCluster(cl)
  result <- t(result)
  ind1 <- which.max(result[,7])
  t3.wei <- matrix(result[ind1,],ncol=1)
  par <- cbind(par,t3.wei)
}
write.csv(par,"par.csv")
quantile(par[1,],c(.025,.975))



##############################
#plot
##############################




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











