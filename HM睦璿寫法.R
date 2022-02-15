# m: burn-in iterations; N:posterior samples; k: thin --------------------------
m=2*10^4;N=2000;k=10

# parameter --------------------------------------------------------------------
a0 <- c(); a1 <- c()
b0 <- c(); b1 <- c()
s0 <- c(); s1 <- c()

# sigma proposal ---------------------------------------------------------------
sigma_proposal=c(10,10,10,10,20,30) 

# acceptance rate --------------------------------------------------------------
alpha1 <- c(); alpha2 <- c(); alpha3 <- c(); alpha4 <- c(); alpha5 <- c(); alpha6 <- c()
update = matrix(NA, nrow = m+N*k, ncol=6)
# the first iterations for using mcmc ------------------------------------------
before <- 10000
stepsize0 <- c(1,1,1,1,1,1)

# initial value ----------------------------------------------------------------
repeat{
  a00 <- runif(1,0,1000); a10 <- runif(1,0,1000)
  b00 <- runif(1,0,1000); b10 <- runif(1,0,1000)
  repeat{
    s00 <- runif(1,0,1000); s10 <- runif(1,-1000,0)
    if(s00+s10*sk[1]>0 & s00+s10*sk[2]>0 & s00+s10*sk[3]>0){break}
  }
  if(alphalikelihood1.1(c(a00*10^(-2),a10*10^(-3),b00*10^(-4),b10*10^(-5),s00,s10),Ti)!= -Inf & 
     is.na(alphalikelihood1.1(c(a00*10^(-2),a10*10^(-3),b00*10^(-4),b10*10^(-5),s00,s10),Ti))==F)
  {break}
}
initial <- c(a00,a10,b00,b10,s00,s10)

start <- Sys.time()
for(i in 1:(m+N*k)){
  # initial value -------------------------------------------------------------------
  a0[1] <- a00; a1[1] <- a10; b0[1] <- b00; b1[1] <- b10; s0[1] <- s00; s1[1] <- s10
  # MH algorithm --------------------------------------------------------------------
  # proposal density -----------------------------------------------------------
  if(i<=before){
    repeat{
      a0[i+1] <- rnorm(1,a0[i],sigma_proposal[1]*stepsize0[1])
      if(a0[i+1]>0){break}
    }
    repeat{
      a1[i+1] <- rnorm(1,a1[i],sigma_proposal[2]*stepsize0[2])
      if(a1[i+1]>0){break}
    }
    repeat{
      b0[i+1] <- rnorm(1, b0[i],sigma_proposal[3]*stepsize0[3])
      if( b0[i+1]>0){break}
    }
    repeat{
      b1[i+1] <- rnorm(1,b1[i],sigma_proposal[4]*stepsize0[4])
      if(b1[i+1]>0 ){break}
    }
    repeat{
      s0[i+1] <- rnorm(1,s0[i],sigma_proposal[5]*stepsize0[5])
      z5.1 <- s0[i+1]+ sk*s1[i]
      if(s0[i+1]>0 & z5.1[1]>0 & z5.1[2]>0 & z5.1[3]>0){break}
    }
    repeat{
      s1[i+1] <- rnorm(1,s1[i],sigma_proposal[6]*stepsize0[6])
      z6.1 <- s0[i]+ sk*s1[i+1]
      if(s1[i+1]<0 & z6.1[1]>0 & z6.1[2]>0 & z6.1[3]>0){break}
    }
    if( i ==2000 | i == 4000 | i == 6000 | i == 8000){
      acpt = c()
      stepchange = c()
      stepsize = c()
      for( k in 1:6){
        acpt[k] = mean(update[1:(i-1),k])                           
        stepchange[k] = log(acpt[k]/(1-acpt[k])) - log(0.4/0.6)           
        stepsize[k] = exp(log(stepsize0[k]) + stepchange[k])
        stepsize0[k] = stepsize[k]
      }
    }
  }else{
    
    repeat{
      a0[i+1] <- rnorm(1,mean(a0[1:i]),sd(a0[1:i]))
      if(a0[i+1]>0){break}
    }
    repeat{
      a1[i+1] <- rnorm(1,mean(a1[1:i]),sd(a1[1:i]))
      if(a1[i+1]>0){break}
    }
    repeat{
      b0[i+1] <- rnorm(1,mean(b0[1:i]),sd(b0[1:i]))
      if(b0[i+1]>0){break}
    }
    repeat{
      b1[i+1] <- rnorm(1,mean(b1[1:i]),sd(b1[1:i]))
      if(b1[i+1]>0 ){break}
    }
    repeat{
      s0[i+1] <- rnorm(1,mean(s0[1:i]),sd(s0[1:i]))
      z5.1 <- s0[i+1]+ sk*s1[i]
      if(s0[i+1]>0 & z5.1[1]>0 & z5.1[2]>0 & z5.1[3]>0){break}
    }
    repeat{
      s1[i+1] <- rnorm(1,mean(s1[1:i]),sd(s1[1:i]))
      z6.1 <- s0[i]+ sk*s1[i+1]
      if(s1[i+1]<0 & z6.1[1]>0 & z6.1[2]>0 & z6.1[3]>0){break}
    }
  }
  # update a0 ------------------------------------------------------------------
  if(i<=before){
    alpha1[i] <- min(alphalikelihood1.1(c(a0[i+1]*10^(-2),a1[i]*10^(-3),b0[i]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti)-
                       alphalikelihood1.1(c(a0[i]*10^(-2),a1[i]*10^(-3),b0[i]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti),0)
  }else{
    alpha1[i] <- min(alphalikelihood1.1(c(a0[i+1]*10^(-2),a1[i]*10^(-3),b0[i]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti)-
                                alphalikelihood1.1(c(a0[i]*10^(-2),a1[i]*10^(-3),b0[i]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti)+
                                dnorm(a0[i+1],mean(a0[1:i]),sd(a0[1:i]),log = T)-dnorm(a0[i],mean(a0[1:i]),sd(a0[1:i]),log = T),0)
  }
  u = runif( 1, 0, 1)
  if(alpha1[i]>log(u)){
    a0[i+1] <- a0[i+1]
    update[i,1] = 1
  }else{
    a0[i+1] <- a0[i]
    update[i,1] = 0
  }
  
  # update a1 ------------------------------------------------------------------
  if(i<=before){
    alpha2[i] <- min(alphalikelihood1.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti)-
                                alphalikelihood1.1(c(a0[i+1]*10^(-2),a1[i]*10^(-3),b0[i]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti),0)
  }else{
    alpha2[i] <- min(alphalikelihood1.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti)-
                                alphalikelihood1.1(c(a0[i+1]*10^(-2),a1[i]*10^(-3),b0[i]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti)+
                                dnorm(a1[i+1],mean(a1[1:i]),sd(a1[1:i]),log = T)-dnorm(a1[i],mean(a1[1:i]),sd(a1[1:i]),log = T),0)
  }
  u = runif( 1, 0, 1)
  if(alpha2[i]>log(u)){
    a1[i+1] <- a1[i+1]
    update[i,2] = 1
  }else{
    a1[i+1] <- a1[i]
    update[i,2] = 0
  }
  
  # update b0 ------------------------------------------------------------------
  if(i<=before){
    alpha3[i] <- min(alphalikelihood2.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti)-
                                alphalikelihood2.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti),0)
  }else{
    alpha3[i] <- min(alphalikelihood2.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti)-
                                alphalikelihood2.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti)+
                                dnorm(b0[i+1],mean(b0[1:i]),sd(b0[1:i]),log = T)-dnorm(b0[i],mean(b0[1:i]),sd(b0[1:i]),log = T),0)
  }
  u = runif( 1, 0, 1)
  if(alpha3[i]>log(u)){
    b0[i+1] <- b0[i+1]
    update[i,3] = 1
  }else{
    b0[i+1] <- b0[i]
    update[i,3] = 0
  }
  
  # update b1 ------------------------------------------------------------------
  if(i<=before){
    alpha4[i] <- min(alphalikelihood2.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i+1]*10^(-5),s0[i],s1[i]),Ti)-
                                alphalikelihood2.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti),0)
  }else{
    alpha4[i] <- min(alphalikelihood2.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i+1]*10^(-5),s0[i],s1[i]),Ti)-
                                alphalikelihood2.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i]*10^(-5),s0[i],s1[i]),Ti)+
                                dnorm(b1[i+1],mean(b1[1:i]),sd(b1[1:i]),log = T)- dnorm(b1[i],mean(b1[1:i]),sd(b1[1:i]),log = T),0)
  }
  u = runif( 1, 0, 1)
  if(alpha4[i]>log(u)){
    b1[i+1] <- b1[i+1]
    update[i,4] = 1
  }else{
    b1[i+1] <- b1[i]
    update[i,4] = 0
  }
  
  # update s0 ------------------------------------------------------------------
  if(i<=before){
    alpha5[i] <- min(alphalikelihood3.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i+1]*10^(-5),s0[i+1],s1[i]),Ti)-
                                alphalikelihood3.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i+1]*10^(-5),s0[i],s1[i]),Ti),0)
  }else{
    alpha5[i] <- min(alphalikelihood3.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i+1]*10^(-5),s0[i+1],s1[i]),Ti)-
                                alphalikelihood3.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i+1]*10^(-5),s0[i],s1[i]),Ti)+
                                dnorm(s0[i+1],mean(s0[1:i]),sd(s0[1:i]),log = T)-dnorm(s0[i],mean(s0[1:i]),sd(s0[1:i]),log = T),0)
  }
  u = runif( 1, 0, 1)
  if(alpha5[i]>log(u)){
    s0[i+1] <- s0[i+1]
    update[i,5] = 1
  }else{
    s0[i+1] <- s0[i]
    update[i,5] = 0
  }
  
  # update s1 ------------------------------------------------------------------
  if(i<=before){
    alpha6[i] <- min(alphalikelihood3.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i+1]*10^(-5),s0[i+1],s1[i+1]),Ti)-
                                alphalikelihood3.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i+1]*10^(-5),s0[i+1],s1[i]),Ti),0)
  }else{
    alpha6[i] <- min(alphalikelihood3.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i+1]*10^(-5),s0[i+1],s1[i+1]),Ti)-
                                alphalikelihood3.1(c(a0[i+1]*10^(-2),a1[i+1]*10^(-3),b0[i+1]*10^(-4),b1[i+1]*10^(-5),s0[i+1],s1[i]),Ti)+
                                dnorm(s1[i+1],mean(s1[1:i]),sd(s1[1:i]),log = T)-dnorm(s1[i],mean(s1[1:i]),sd(s1[1:i]),log = T),0)
  }
  u = runif( 1, 0, 1)
  if(alpha6[i]>log(u)){
    s1[i+1] <- s1[i+1]
    update[i,6] = 1
  }else{
    s1[i+1] <- s1[i]
    update[i,6] = 0
  }# end -----------------------------------------------------------------------
}
end <- Sys.time()
end-start