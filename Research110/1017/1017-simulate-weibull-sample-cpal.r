#-------------- calculate coverage probability & average length -----------------#

# Above is one sample of simulation
# We need to simulate large sample

#-----------------------------------------------------------------------------------

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

f <- function(g){
  #likelihood fn 
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
    ans= g*lgamma(1/c+1)+g*log(c)+g*log(a)+sum(k)-sum(k1)+b*sumt
    -ans
  }
  
  ll1.wei <- function(theta,data){
    ll.wei(theta[1]*10^(-1),theta[2]*10^(-1),theta[3]*10^(-1),data)
  }
  
  x <- apply(matrix(est0[,3],ncol=1),1,
             function(x){rweibull(g,shape=x,scale=(1/gamma(1/x+1)))})
  x1 <- apply(x,2,function(s){cumsum(s)})
  x2 <- matrix(0,ncol=3)
  for(i in 1:nk){
    t <- cbind (rep(est0[i,1],g),rep(est0[i,2],g),x1[,i])
    x2 <- rbind(x2,t)
  }
  x2 <- x2[-1,]
  Ti<- apply(x2,1,function(x){log(x[2]*x[3]/x[1]+1)/x[2]})
  Ti <- matrix(Ti,ncol=nk)
  Ti <- rbind(rep(0,nk),Ti)
  rm(i,t)
  #estimate
  r1 <- runif(10,0,100);r2 <- runif(10,0,100);r3 <- runif(10,0,100)
  table2 <- matrix(0,nrow = 9,ncol =7)
  for(k in 1:9){
    optim.result=matrix(0,nrow = 1,ncol = 7)
    for (i in 1:10){
      fit=try(optim(c(r1[i],r2[i],r3[i]),ll1.wei,hessian=T,data=Ti[,k]))
      if(class(fit)=="try-error"){
        next
      }
      fit1=try(optim(fit$par,ll1.wei,hessian=T,data=Ti[,k]))
      if(class(fit1)=="try-error"){
        next
      }
      if (abs(fit1$value-fit$value)>10^(-3)){
        fit <- fit1
        fit1<- optim(fit$par,ll1.wei,hessian=T,data=Ti[,k])
      }
      I=solve(fit1$hessian)
      sda <- sqrt(I[1,1]*10^(-2))
      sdb <- sqrt(I[2,2]*10^(-2))
      sds <- sqrt(I[3,3]*10^(-2))
      t=cbind(fit1$par[1]*10^(-1),sda,
              fit1$par[2]*10^(-1),sdb,
              abs(fit1$par[3]*10^(-1)),sds,
              (-fit1$value))
      optim.result=rbind(optim.result,t)
    }
    optim.result <- optim.result[-1,]
    ind <- which.max(optim.result[,7])
    table2[k,] <- optim.result[ind,]
  }
  ci <- cbind(t(apply(table2[,1:2],1,function(x){x[1]+c(-1,1)*1.96*x[2]})),
              t(apply(table2[,3:4],1,function(x){x[1]+c(-1,1)*1.96*x[2]})),
              t(apply(table2[,5:6],1,function(x){x[1]+c(-1,1)*1.96*x[2]})))
  cp <- matrix(0,ncol =3,nrow = 9)
  al <- matrix(0,ncol =3,nrow = 9)
  
  for(j in 1:9){
    for(h in 1:3){
      theoretical <- est0[j,h]
      left <- ci[j,2*h-1]
      right <- ci[j,2*h]
      if(left < theoretical  & theoretical  < right){
        cp[j,h] <- 1;al[j,h] <- right-left
      }else{
        cp[j,h] <- 0;al[j,h] <- right-left
      }
    }
  }
  cp <- as.numeric(cp)
  al <- as.numeric(al)
  return(c(cp,al))
}




start_time <- Sys.time()
result45 <- replicate(1000,f(45))
end_time <- Sys.time()
end_time - start_time

cpal45 <- matrix(apply(result45,1,function(x){mean(x)}),ncol=6)
cpal45 <- cbind(cpal45[,1],cpal45[,4],cpal45[,2],cpal45[,5],cpal45[,3],cpal45[,6])
colnames(cpal45) <- c("a.CP","AL","b.CP","AL","s.CP","AL")
write.csv(cpal45,"cpal45.csv")



result100 <- c()
start_time <- Sys.time()
for (i in 1:500){
  t=try(f(100))
  if(class(t)=="try-error"){
    next
  }
  result100 <- append(result100,t)
}
end_time <- Sys.time()
end_time - start_time

result100 <- matrix(result100,nrow = 48)

ncol(result100)
write.csv(result100,"result100.csv")

cpal100 <- matrix(apply(result100,1,function(x){mean(x)}),ncol=6)
cpal100 <- cbind(cpal100[,1],cpal100[,4],cpal100[,2],cpal100[,5],cpal100[,3],cpal100[,6])
colnames(cpal100) <- c("a.CP","AL","b.CP","AL","s.CP","AL")
write.csv(cpal100,"cpal100.csv")














