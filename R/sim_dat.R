library(gtools)
library(simsurv)

sim_dat <- function(M,p,ni,k,ksurv,beta){
  
  W <- matrix(rexp(p*k,2),nrow=p,ncol=k)
  
  for(i in 1:k){
    W[((i-1)*25+1):(i*25),i] <- 10
  }
  
  H <- list()
  for(i in 1:M){
    H[[i]] <- matrix(rexp(ni*k,1/.25),nrow=k,ncol=ni)
  }
  Hbind <- do.call('cbind',H)
  
  # si <- sim.survdata(N=ni, T=80, num.data.frames = 1, fixed.hazard = TRUE,
  #              X=t(Hbind[1:ksurv,]),beta=beta)
  xsurv <- data.frame(t(Hbind)[,1:ksurv])
  names(beta) <- colnames(xsurv)
  si <- simsurv(dist = 'weibull',lambdas = .0001, gammas = 2.1, x = xsurv, betas = beta, maxt=60)
  
  y <- list()
  delta <- list()
  X <- list()
  for(i in 1:M){
    samps <- si[((i-1)*ni+1):(ni*i),]
    y[[i]] <- samps$eventtime
    delta[[i]] <- samps$status
    X[[i]] <- W %*% H[[i]]
  }
  
  return(list(X=X,y=y,delta=delta,Wtrue=W,Htrue=H))
}





