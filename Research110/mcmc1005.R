
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
tll <- function(theta,data){ # ll is -loglikeihood, so here is also -logL
  a0 <- theta[1];a1 <- theta[2];b0 <- theta[3];b1 <- theta[4];s0 <- theta[5];s1 <- theta[6]
  aS <- apply(matrix(c(0,1),ncol = 1),1,function(x) a0+a1*x)
  bS <- apply(matrix(c(0,1),ncol = 1),1,function(x) b0+b1*x)
  sS <- apply(matrix(c(0,1),ncol = 1),1,function(x) s0+s1*x)
  t1 <- sum(apply(matrix(1,ncol = 1),1,function(x) ll(aS[1],bS[1],sS[1],data[,x])))
  t2 <- sum(apply(matrix(2,ncol = 1),1,function(x) ll(aS[2],bS[2],sS[2],data[,x])))
  ans <- t1+t2
  -ans
}

m=10^4;N=2000;k=10
a00 <- 0.1;a10 <- 0.1;b00 <- 0.1;b10 <- 0.1;s00 <-0.1;s10 <- 0.1;sigma=0.01
a0 <- c();a1 <- c();b0 <- c();b1 <- c();s0 <- c();s1 <- c()
u <- runif(m+N*k)
start <- Sys.time()
for(i in 1:(m+N*k)){
  repeat{
    z1 <- rnorm(1,a00,sigma)
    if(z1>0 & z1<1){break}
  }
  p1 <- exp(-(z1-0.0916)^2/(2*sigma^2))
  p0 <- exp(-(a00-0.0916)^2/(2*sigma^2))
  alpha1 <- min(exp(tll(c(z1,a10,b00,b10,s00,s10),Ti)-tll(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha1>u[i]){
    a0[i] <- z1
  }else{
    a0[i] <- a00
  }
  a00 <- a0[i]
  repeat{
    z2 <- rnorm(1,a10,sigma)
    if(z2>0 & z2<1){break}
  }
  p1 <- exp(-(z2-0.000836)^2/(2*sigma^2))
  p0 <- exp(-(a10-0.000836)^2/(2*sigma^2))
  alpha2 <- min(exp(tll(c(a00,z2,b00,b10,s00,s10),Ti)-tll(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha2>u[i]){
    a1[i] <- z2
  }else{
    a1[i] <- a10
  }
  a10 <- a1[i]
  repeat{
    z3 <- rnorm(1,b00,sigma)
    if(z3>0 & z3<1){break}
  }
  p1 <- exp(-(z3-0.000203)^2/(2*sigma^2))
  p0 <- exp(-(b00-0.000203)^2/(2*sigma^2))
  alpha3 <- min(exp(tll(c(a00,a10,z3,b10,s00,s10),Ti)-tll(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha3>u[i]){
    b0[i] <- z3
  }else{
    b0[i] <- b00
  }
  b00 <- b0[i]
  repeat{
    z4 <- rnorm(1,b10,sigma)
    if(z4>0 & z4<1){break}
  }
  p1 <- exp(-(z4-0.0000165)^2/(2*sigma^2))
  p0 <- exp(-(b10-0.0000165)^2/(2*sigma^2))
  alpha4 <- min(exp(tll(c(a00,a10,b00,z4,s00,s10),Ti)-tll(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha4>u[i]){
    b1[i] <- z4
  }else{
    b1[i] <- b10
  }
  b10 <- b1[i]
  repeat{
    z5 <- rnorm(1,s00,sigma)
    if(z5>0 & z5<1){break}
  }
  p1 <- exp(-(z5-0.00291)^2/(2*sigma^2))
  p0 <- exp(-(s00-0.00291)^2/(2*sigma^2))
  alpha5 <- min(exp(tll(c(a00,a10,b00,b10,z5,s10),Ti)-tll(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha5>u[i]){
    s0[i] <- z5
  }else{
    s0[i] <- s00
  }
  s00 <- s0[i]
  repeat{
    z6 <- rnorm(1,s10,sigma)
    if(z6>0 & z6<1){break}
  }
  p1 <- exp(-(z6-0.00599)^2/(2*sigma^2))
  p0 <- exp(-(s10-0.00599)^2/(2*sigma^2))
  alpha6 <- min(exp(tll(c(a00,a10,b00,b10,s00,z6),Ti)-tll(c(a00,a10,b00,b10,s00,s10),Ti))*p1/p0,1)
  if(alpha6>u[i]){
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



