posterior <- function(theta,y){(theta^(y-0.5))*((1-theta)^(50-y-0.5))}
f <- function(m=1000,N=1000,k=30,theta0=0.5,sigma){
  a <-c();theta <- c()
  y0 <- rbinom(1,50,0.8)
  u <- runif(m+N*k)
  for(i in 1:(m+N*k)){
    repeat{
      z <- rnorm(1,theta0,sigma)
      if(z>0 & z<1){break}
    }
    a[i] <- alpha <- min(posterior(z,y0)/posterior(theta0,y0),1)
    if(alpha>u[i]){
      theta[i] <- z
    }else{
      theta[i] <- theta0
    }
    theta0 <- theta[i]
  }
  return(list(a,u,theta))
}
set.seed(109225017)
result1 <- f(sigma=0.1)
m=1000;N=1000;k=30
result1.theta <- result1[[3]][seq(m+k,m+N*k,k)]
table <- matrix(c(mean(result1.theta),var(result1.theta),sum(result1[[1]]>result1[[2]])/(m+N*k)),nrow=1)
set.seed(109225017)
result2 <- f(sigma=0.5)
result2.theta <- result2[[3]][seq(m+k,m+N*k,k)]
table <- rbind(table,c(mean(result2.theta),var(result2.theta),sum(result2[[1]]>result2[[2]])/(m+N*k)))
colnames(table) <- c("E(theta|y)","Var(theta|y)","AR")
rownames(table) <- c("sigma=0.1","sigma=0.5")
table <- as.data.frame(table)

par(mfrow=c(3,2))
plot(1:(m+N*k),result1[[3]],type="l",xlab = "iteration number",ylab="MCMC sample",main=expression
     (sigma==0.1))
abline(v=m,col="red",lty=2)

plot(1:(m+N*k),result2[[3]],type="l",xlab = "iteration number",ylab="MCMC sample",main=expression
     (sigma==0.5))
abline(v=m,col="red",lty=2)

plot(1:N,result1.theta,type="l",xlab = "iteration number",ylab="MCMC sample",main=expression(sigma==0.1))
abline(h=0.8,col="red",lty=2)

plot(1:N,result2.theta,type="l",xlab = "iteration number",ylab="MCMC sample",main=expression(sigma==0.5))
abline(h=0.8,col="red",lty=2)

hist(result1.theta,freq = F,xlim=c(0,1),xlab=expression(theta),main=expression(sigma==0.1))
curve(dbeta(x,0.5,0.5),add=T,col="red")

hist(result2.theta,freq = F,xlim=c(0,1),xlab=expression(theta),main=expression(sigma==0.5))
curve(dbeta(x,0.5,0.5),add=T,col="red")


#2
y0 <- 0.5
N = m = 1000;k <- 30
x <-c(); y <- c()
set.seed(109225017)
for(i in 1:(m+N*k)){
  x[i] <- rbinom(1,20,y0)
  y[i] <- rbeta(1,x[i]+2,20-x[i]+4)
  y0 <- y[i]
}
result <- y[seq(m+k,m+N*k,k)]
mean(result);var(result)
xaxis <- seq(0.01,0.99,0.02)
fyhat <- sapply(xaxis, function(y){
  mean(dbeta(y,x[seq(m+k,m+N*k,k)]+2,20-x[seq(m+k,m+N*k,k)]+4))
})
plot(xaxis,fyhat,type = "l",xlab = "y",ylab =expression(hat(f)))
