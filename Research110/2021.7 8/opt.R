Opt<-function(x){
  iter=x[7]
  initial=c(x[1],...,x[6])
  reeor=1
  before=0
  while(abs(error)-1e-15>0){
    step1=try(optim(initial,likelihood,control=list(factr=1e-15,maxit=10^9)),silent = T)
    if(step1$convergence==0){
      initial <- step1$par
      error <- step1$value-before
      before <- step1$value
    }else{
      step1$value <- Inf
      return(step1
             )
    }
    print(c(iter,-before,step1$convergence))
  }
  return(step1)
  
}
apply(as.matrix(rn),1,function(x){tryCatch(Opt(x),error=function (e) NA)})