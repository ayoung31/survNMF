library(NMF)
init_H <- function(X,k,method='IndNMF',ninit=50){
  
  if(method=='IndNMF'){
    M <- length(X)
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
  }
  return(H0)
}

H0 <- init_H(X,k)
