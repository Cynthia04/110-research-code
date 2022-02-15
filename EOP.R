
#### EOP #######
#已標準化應力

#-------------------------------------- 
# renewal distribution F = N(1,sigma^2) 
#-------------------------------------- 

# # No standardize stress
# t3est <- c(9.14*10^(-2),1.02*10^(-3),2.1*10^(-4),1.12*10^(-5),2.92*10^(-3),4.42*10^(-4))
# t3reg <- c(9.15*10^(-2),1*10^(-3),2.07*10^(-4),1.2*10^(-5),1.72*10^(-3),3.46*10^(-4))
# est <- t3reg 
# S0 <- c()
# for (i in c(1,3,5)){
#   t <- est[i]+0.5*est[i+1]
#   S0 <- append(S0,t)
# }
# EZ <- function(i){
#   a <- S0[1];b <- S0[2];s <- S0[3]
#   k <- a/b
#   1/b*(log((i+k)/(i-1+k))+s^2/2*((i-1)/(i-1+k)^2-i/(i+k)^2))-8
# }
# ceiling(uniroot(EZ, lower = 1, upper = 200)$root)

# standardize stress
t3.nor <- read.csv("t3_nor.csv")
theta <- as.numeric(t3.nor[1:6,1])
# EZ <- function(x,i){
#   a <- x[1];b <- x[2];c <- exp(x[3])
#   mu_Y <- function(i) {1+i*b/a}
#   Var_Y <- function(i) {(b/a)^2*i*(exp(c^2)-1)}
#   1/b*((log(mu_Y(i))-Var_Y(i)/(2*mu_Y(i)^2))-(log(mu_Y(i-1))-Var_Y(i-1)/(2*mu_Y(i-1)^2)))-8
# }
EZ <- function(x,i){
  a <- x[1];b <- x[2];s <- x[3]
  k <- a/b
  1/b*(log((i+k)/(i-1+k))+s^2/2*((i-1)/(i-1+k)^2-i/(i+k)^2))-8
}
EZ.nor <- function(i){
  EZ(theta[c(1,3,5)],i)
}
ceiling(uniroot(EZ.nor, lower = 1, upper = 200)$root)#154

#------------------------------ sigma=1 ---------------------------------------#
t3.nor_sigma <- read.csv("t3_nor_sigma.csv")
theta <- as.numeric(t3.nor_sigma[1:6,1])
EZ <- function(x,i){
  a <- x[1];b <- x[2];c <- exp(x[3]) ### !!!
  mu_Y <- function(i) {1+i*b/a*c}
  Var_Y <- function(i) {(b/a)^2*i}
  1/b*((log(mu_Y(i))-Var_Y(i)/(2*mu_Y(i)^2))-(log(mu_Y(i-1))-Var_Y(i-1)/(2*mu_Y(i-1)^2)))-8
}
EZ.nor_sigma <- function(i){
  EZ(theta[c(1,3,5)],i)
}
ceiling(uniroot(EZ.nor_sigma, lower = 1, upper = 200)$root) #155

#-------------------------------------------------------  
# renewal distribution F = Lognormal(-alpha^2/2,alpha^2) 
#------------------------------------------------------- 
t3.log <- read.csv("t3_log.csv")
theta <- as.numeric(t3.log[1:6,1])
EZ <- function(x,i){
  a <- x[1];b <- x[2];c <- x[3]
  mu_Y <- function(i) {1+b*i/a}
  Var_Y <- function(i) {(b/a)^2*(exp(c^2)-1)*i}
  1/b*((log(mu_Y(i))-Var_Y(i)/(2*mu_Y(i)^2))-(log(mu_Y(i-1))-Var_Y(i-1)/(2*mu_Y(i-1)^2)))-8
}
EZ.log <- function(i){
  EZ(theta[c(1,3,5)],i)
}
ceiling(uniroot(EZ.log, lower = 1, upper = 200)$root)#154

#--------------------------------- mu=1 ---------------------------------------#
t3.log_mu <- read.csv("t3_log_mu.csv")
theta <- as.numeric(t3.log_mu[1:6,1])
EZ <- function(x,i){
  a <- x[1];b <- x[2];c <- x[3]
  mu_Y <- function(i) {1+i*b/a*exp(1+c^2/2)}
  Var_Y <- function(i) {(b/a)^2*((exp(c^2)-1)*exp(2+c^2))*i}
  1/b*((log(mu_Y(i))-Var_Y(i)/(2*mu_Y(i)^2))-(log(mu_Y(i-1))-Var_Y(i-1)/(2*mu_Y(i-1)^2)))-8
}
EZ.log_mu <- function(i){
  EZ(theta[c(1,3,5)],i)
}
ceiling(uniroot(EZ.log_mu, lower = 1, upper = 200)$root) #154


#--------------------------------------------------------- 
# renewal distribution F = Weibull(1/gamma(1+1/beta),beta) 
#--------------------------------------------------------- 
t3.wei <- read.csv("mle.csv")[,-1]
theta <- as.numeric(t3.wei[1:6])
# EZ <- function(x,i){
#   a <- x[1];b <- x[2];c <- x[3]
#   mu_Y <- function(i) {1+b*i/a}
#   Var_Y <- function(i) {(b/a)^2*((1/gamma(1/c+1))^2*gamma(2/c+1)-1)*i}
#   1/b*((log(mu_Y(i))-Var_Y(i)/(2*mu_Y(i)^2))-(log(mu_Y(i-1))-Var_Y(i-1)/(2*mu_Y(i-1)^2)))-8
# }
EZ <- function(x,i){
  a <- x[1];b <- x[2];c <- x[3]
  k <- a/b; var <- (1/gamma(1/c+1))^2*gamma(2/c+1)-1
  1/b*(log((i+k)/(i-1+k))+var/2*((i^2-i-k^2)/((i-1+k)^2*(i+k)^2)))-8
}
EZ.wei <- function(i){
  EZ(theta[c(1,3,5)],i)
}
ceiling(uniroot(EZ.wei, lower = 1, upper = 200)$root)#156

#------------------------------ eta=1 -----------------------------------------#
t3.wei_eta <- read.csv("t3.wei_eta_1.csv")[,-1]
theta <- as.numeric(t3.wei_eta[1:6])
EZ <- function(x,i){
  a <- x[1];b <- x[2];c <- x[3]
  mu_Y <- function(i) {1+b*i/a*gamma(1/c+1)}
  Var_Y <- function(i) {(b/a)^2*(gamma(2/c+1)-gamma(1/c+1)^2)*i}
  1/b*((log(mu_Y(i))-Var_Y(i)/(2*mu_Y(i)^2))-(log(mu_Y(i-1))-Var_Y(i-1)/(2*mu_Y(i-1)^2)))-8
}
EZ.wei_eta <- function(i){
  EZ(theta[c(1,3,5)],i)
}
ceiling(uniroot(EZ.wei_eta, lower = 1, upper = 200)$root)#156






