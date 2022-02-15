# ---------------------------- Trace of mcmc sample ----------------------------
traceplot <- function(data,ylab0,main0){
  plot(data,type = "l",xlab = "iteration times (40000)",
       ylab=ylab0,main=main0,
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
}
par(mfrow=c(3,2))
traceplot(a0,expression(a[0]),expression("Trace of MCMC sample of "~ a[0]))
traceplot(a1,expression(a[1]),expression("Trace of MCMC sample of "~ a[1]))
traceplot(b0,expression(b[0]),expression("Trace of MCMC sample of "~ b[0]))
traceplot(b1,expression(b[1]),expression("Trace of MCMC sample of "~ b[1]))
traceplot(s0,expression(beta[0]),expression("Trace of MCMC sample of "~ beta[0]))
traceplot(s1,expression(beta[1]),expression("Trace of MCMC sample of "~ beta[1]))
par(mfrow=c(1,1))
# ----------------------- Histogram of posterior sample ------------------------
posterior.dist <- function(data,bre,adj,abv,xname,main0){
  hist(data,breaks = bre,
       col="gray90", # column color
       border="black",
       prob = TRUE, # show densities instead of frequencies
       xlab = xname,
       main=main0,
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  lines(density(data, adjust = adj),
        lwd = 2, # thickness of line
        col = "red")
  abline(v= abv, lwd = 2, lty=2, col = "blue")
}
par(mfrow=c(3,2))
f <- function(x){
  posterior.dist(x[,1],8,5,t3.wei[1,1],expression(a[0]),
                 expression("Histogram of posterior sample of "~a[0]))
  posterior.dist(x[,2],10,3,t3.wei[2,1],expression(a[1]),
                 expression("Histogram of posterior sample of "~a[1]))
  posterior.dist(x[,3],10,3,t3.wei[3,1],expression(b[0]),
                 expression("Histogram of posterior sample of "~b[0]))
  posterior.dist(x[,4],10,3,t3.wei[4,1],expression(b[1]),
                 expression("Histogram of posterior sample of "~b[1]))
  posterior.dist(x[,5],10,3,t3.wei[5,1],expression(beta[0]),
                 expression("Histogram of posterior sample of "~beta[0]))
  posterior.dist(x[,6],10,3,t3.wei[6,1],expression(beta[1]),
                 expression("Histogram of posterior sample of "~beta[1]))
}
par(mfrow=c(3,2))
f(posterior_sample1)
par(mfrow=c(1,1))

par(mfrow=c(3,2))
f(posterior_sample2)
par(mfrow=c(1,1))

par(mfrow=c(3,2))
f(posterior_sample3)
par(mfrow=c(1,1))

par(mfrow=c(3,2))
f(posterior_sample)
par(mfrow=c(1,1))
