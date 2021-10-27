#--------------- Use the normal estimate of table 2 to simulate data -----------#
est <- data.frame(a=c(9.219,9.269,9.235,9.422,9.503,9.434,9.649,9.638),
                  b=c(2.182,2.215,2.258,2.435,2.394,2.449,2.771,2.583),
                  sigma=c(2.339,2.080,1.782,3.357,2.518,2.541,3.944,4.053))
est1 <- cbind(est[,1]*10^(-2),est[,2]*10^(-4),est[,3]*10^(-3))

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
  
  ll1 <- function(theta,data){
    ll(theta[1]*10^(-3),theta[2]*10^(-3),theta[3]*10^(-3),data)
  }
  
  x <- apply(matrix(c(est[,3]*10^(-3)),ncol = 1),1,function(x){rnorm(g,mean=1,sd=x)})
  x1 <- apply(x,2,function(s){cumsum(s)})
  x2 <- matrix(0,ncol=3)
  for(i in 1:8){
    t <- cbind (rep(est[i,1],g),rep(est[i,2],g),x1[,i])
    x2 <- rbind(x2,t)
  }
  x2 <- x2[-1,]
  Ti<- apply(x2,1,function(x){log(x[2]*x[3]*10^(-2)/x[1]+1)/(x[2]*10^(-4))})
  Ti <- matrix(Ti,ncol=8)
  Ti <- rbind(rep(0,8),Ti)
  #estimate
  r1 <- runif(10,0,100);r2 <- runif(10,0,100);r3 <- runif(10,0,100)
  table2 <- matrix(0,nrow = 8,ncol =7)
  for(k in 1:8){
    optim.result=matrix(0,nrow = 1,ncol = 7)
    for (i in 1:10){
      fit=try(optim(c(r1[i],r2[i],r3[i]),ll1,hessian=T,data=Ti[,k]))
      if(class(fit)=="try-error"){
        next
      }
      fit1=try(optim(fit$par,ll1,hessian=T,data=Ti[,k]))
      if(class(fit1)=="try-error"){
        next
      }
      if (abs(fit1$value-fit$value)>10^(-3)){
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
    ind <- which.max(optim.result[,7])
    table2[k,] <- optim.result[ind,]
  }
  ci <- cbind(t(apply(table2[,1:2],1,function(x){x[1]+c(-1,1)*1.96*x[2]})),
              t(apply(table2[,3:4],1,function(x){x[1]+c(-1,1)*1.96*x[2]})),
              t(apply(table2[,5:6],1,function(x){x[1]+c(-1,1)*1.96*x[2]})))
  cp <- matrix(0,ncol =3,nrow = 8)
  al <- matrix(0,ncol =3,nrow = 8)

  for(j in 1:8){
    for(h in 1:3){
      theoretical <- est1[j,h]
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
result45 <- replicate(500,f(45))
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












