library(gtools)
library(simsurv)

sim_dat <- function(M,p,ni,k,lb,lf,ksurv,beta){
  
  W <- matrix(rexp(p*k,lb),nrow=p,ncol=k)
  
  genes <- list()
  for(i in 1:k){
    genes[[i]] <- sample(1:nrow(W),25,replace=FALSE)
    W[genes[[i]],i] <- rexp(25,lf)
    
  }
  
  #set first ksurv factors to small values compared to others
  H <- list()
  
  for(i in 1:M){
    H[[i]] <- matrix(nrow=k,ncol=ni)
    H[[i]][1:ksurv,] <- rexp(ni*ksurv,1/.08)
  }
  
  for(i in 1:M){
    H[[i]][(ksurv+1):k,] <- matrix(rexp(ni*(k-ksurv),1/.25),nrow=(k-ksurv),ncol=ni)
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
  
  return(list(X=X,y=y,delta=delta,Wtrue=W,Htrue=H,gTrue=genes))
}





