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

m=10^4;N=2000;k=10
a00 <- 3;a10 <- 2;b00 <- 5;b10 <- 4;s00 <-20;s10 <- -5;sigma=0.1
a0 <- c();a1 <- c();b0 <- c();b1 <- c();s0 <- c();s1 <- c()
alpha1 <- c(); alpha2 <- c(); alpha3 <- c(); alpha4 <- c(); alpha5 <- c(); alpha6 <- c()
u <- runif(m+N*k)
start <- Sys.time()
for(i in 1:(m+N*k)){
  #a0
  z1 <- rnorm(1,a00,sigma)
  p1 <- exp(-(z1-3.00162)^2/(2*sigma^2))
  p0 <- exp(-(a00-3.00162)^2/(2*sigma^2))
  alpha1[i] <- alpha <- min(exp(tll.wei(c(z1,a10,b00,b10,s00,s10),Ti)-tll.wei(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha>u[i]){
    a0[i] <- z1
  }else{
    a0[i] <- a00
  }
  a00 <- a0[i]
  #a1
  z2 <- rnorm(1,a10,sigma)
  p1 <- exp(-(z2-2.00355)^2/(2*sigma^2))
  p0 <- exp(-(a10-2.00355)^2/(2*sigma^2))
  alpha2[i] <- alpha <- min(exp(tll.wei(c(a00,z2,b00,b10,s00,s10),Ti)-tll.wei(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha>u[i]){
    a1[i] <- z2
  }else{
    a1[i] <- a10
  }
  a10 <- a1[i]
  #b0
  z3 <- rnorm(1,b00,sigma)
  p1 <- -(z3-4.99938)^2/(2*sigma^2)
  p0 <- -(b00-4.99938)^2/(2*sigma^2)
  alpha3[i] <- alpha <- min(exp(tll.wei(c(a00,a10,z3,b10,s00,s10),Ti)-tll.wei(c(a00,a10,b00,b10,s00,s10),Ti))*exp(p1-p0),1)
  if(alpha>u[i]){
    b0[i] <- z3
  }else{
    b0[i] <- b00
  }
  b00 <- b0[i]
  #b1
  z4 <- rnorm(1,b10,sigma)
  p1 <- exp(-(z4-3.9979)^2/(2*sigma^2))
  p0 <- exp(-(b10-3.9979)^2/(2*sigma^2))
  alpha4[i] <- alpha <- min(exp(tll.wei(c(a00,a10,b00,z4,s00,s10),Ti)-tll.wei(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha>u[i]){
    b1[i] <- z4
  }else{
    b1[i] <- b10
  }
  b10 <- b1[i]
  #beta0
  z5 <- rnorm(1,s00,sigma)
  p1 <- -(z5-20.23702)^2/(2*sigma^2)
  p0 <- -(s00-20.23702)^2/(2*sigma^2)
  alpha5[i] <- alpha <- min(exp(tll.wei(c(a00,a10,b00,b10,z5,s10),Ti)-tll.wei(c(a00,a10,b00,b10,s00,s10),Ti))*exp(p1-p0),1)
  if(alpha>u[i]){
    s0[i] <- z5
  }else{
    s0[i] <- s00
  }
  s00 <- s0[i]
  #beta1
  z6 <- rnorm(1,s10,sigma)
  p1 <- -(z6+5.08808)^2/(2*sigma^2)
  p0 <- -(s10+5.08808)^2/(2*sigma^2)
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

mcmc <- function(x){return(c(mean(x),sd(x),quantile(x, probs = c(0.025, 0.975))))}
mcmc_table <- t(apply(rbind(a0[seq(m+k,m+N*k,k)],a1[seq(m+k,m+N*k,k)],b0[seq(m+k,m+N*k,k)],
                            b1[seq(m+k,m+N*k,k)],s0[seq(m+k,m+N*k,k)],s1[seq(m+k,m+N*k,k)]),1,mcmc))

mcmc_table <- cbind(c(3,2,5,4,20,-5),mcmc_table)
colnames(mcmc_table) <- c("True","mean","sd","2.5%","97.5%")
mcmc_table

par(mfrow=c(3,2))
plot(a0[seq(k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of a0")
plot(a1[seq(k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of a1")
plot(b0[seq(k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of b0")
plot(b1[seq(k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of b1")
plot(s0[seq(k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of sigma0")
plot(s1[seq(k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of sigma1")
par(mfrow=c(1,1))

length(seq(m+k,m+N*k,k))

par(mfrow=c(3,2))
plot(a0[seq(m+k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of a0")
plot(a1[seq(m+k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of a1")
plot(b0[seq(m+k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of b0")
plot(b1[seq(m+k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of b1")
plot(s0[seq(m+k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of sigma0")
plot(s1[seq(m+k,m+N*k,k)],type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of sigma1")
par(mfrow=c(1,1))

mcmc_table
t3
write.csv(mcmc_table,"mcmc.csv")




par(mfrow=c(3,2))
plot(a0,type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of a0")
plot(a1,type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of a1")
plot(b0,type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of b0")
plot(b1,type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of b1")
plot(s0,type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of sigma0")
plot(s1,type = "p",xlab = "iteration number",ylab="MCMC sample",main="posterior sample of sigma1")
par(mfrow=c(1,1))



