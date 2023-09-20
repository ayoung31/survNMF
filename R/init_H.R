library(NMF)

#' @export
init_H <- function(X,k,method='merge',ninit=10){
  
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
  }else if(method=='merge'){
    fit <- nmf(do.call('cbind',X),k,nrun=ninit)
    H <- fit@fit@H %*% diag(1/colSums(fit@fit@H))
    H0 <- list()
    nc <- ncol(X[[1]])
    H0[[1]] <- H[,1:nc]
    for(i in 2:M){
      ncprev <- ncol(X[[i-1]])
      nc <- ncol(X[[i]])
      H0[[i]] <- H[,(ncprev+1):(ncprev+nc)]
    }
  }
  return(H0)
}
