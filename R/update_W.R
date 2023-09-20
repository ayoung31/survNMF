#' @export
update_W <- function(X,H,theta){
  M <- length(X)
  a <- list()
  b <- list()
  for(i in 1:M){
    a[[i]] <- sqrt(theta[i]) * H[[i]] %*% t(H[[i]])
    b[[i]] <- sqrt(theta[i]) * X[[i]] %*% t(H[[i]])
  }
  A <- t(Reduce('+',a))
  B <- t(Reduce('+',b))
  
  fit <- .fcnnls(A,B)
  return(fit$coef)
}