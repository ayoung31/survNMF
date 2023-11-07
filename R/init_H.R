
#' @export
init_H <- function(X,k,y,delta,theta,alpha,lambda,eta,method='random',ninit=5){
  M <- length(X)
  if(method=='IndNMF'){
    
    H0 <- list()
    for(i in 1:M){
      fit <- list()
      resid <- numeric()
      for(j in 1:ninit){
        fit[[j]] <- nmf(X[[i]],k)
        resid[j] <- fit[[j]]@residuals
      }
      which.min(resid)
      best.fit <- fit[[which.min(resid)]]
      H0[[i]] <- best.fit@fit@H %*% diag(1/colSums(best.fit@fit@H))
    }
  }else if(method=='merge'){
    fit <- NMF::nmf(do.call('cbind',X),k,nrun=ninit)
    H <- fit@fit@H %*% diag(1/colSums(fit@fit@H))
    H0 <- list()
    nc <- ncol(X[[1]])
    H0[[1]] <- H[,1:nc]
    if(M>1){
      for(i in 2:M){
        ncprev <- ncol(X[[i-1]])
        nc <- ncol(X[[i]])
        H0[[i]] <- H[,(ncprev+1):(ncprev+nc)]
      }
    }
  }else if(method=='random'){
    cprev <- 0
    for(j in 1:ninit){
      set.seed(j)
      H0 <- list()
      for(i in 1:M){
        n <- ncol(X[[i]])
        H0[[i]] <- matrix(runif(k*n),nrow=k,ncol=n)
      }
      #call optimize loss
      fit <- optimize_loss(X,H0,k,y,delta,theta,alpha,lambda,eta,maxit=10)
      Hcurr <- t(do.call('cbind',fit$H))
      cind <- cvwrapr::getCindex(Hcurr %*% fit$beta,Surv(unlist(y),unlist(delta)))
      if(cind > cprev){
        cprev <- cind
        best_fit <- fit
      }
    }
    H0 <- best_fit$H
  }
  
  
  return(H0)
}
