library(gtools)
library(coxed)

sim_data <- function(M,p,ni,k,ksurv,beta){
  
  W <- matrix(.0001,nrow=p,ncol=k)
  
  for(i in 1:k){
    W[((i-1)*25+1):(i*25),i] <- 10
  }
  
  H <- list()
  for(i in 1:M){
    H[[i]] <- t(rdirichlet(ni,rep(1,k)))
  }
  Hbind <- do.call('cbind',H)
  
  si <- sim.survdata(N=ni, T=80, num.data.frames = 1, fixed.hazard = TRUE,
               X=t(Hbind[1:ksurv,]),beta=beta)
  
  y <- list()
  delta <- list()
  X <- list()
  for(i in 1:M){
    samps <- si$data[((i-1)*ni+1):(ni*i),]
    y[[i]] <- samps$y
    delta[[i]] <- as.numeric(samps$failed)
    X[[i]] <- W %*% H[[i]]
  }
  
  return(list(X=X,y=y,delta=delta,Wtrue=W,Htrue=H))
}





