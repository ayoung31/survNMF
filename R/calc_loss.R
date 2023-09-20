#' @export
calc_loss <- function(X,W,H,beta,alpha,y,delta,theta,lambda){
  M <- length(X)
  p <- 0
  for(i in 1:M){
    ni <- ncol(X[[i]])
    yi <- y[[i]]
    Xi <- X[[i]]
    Hi <- H[[i]]
    deltai <- delta[[i]]
    p1 <- norm(Xi-W%*%Hi,'F')^2
    p2 <- 0
    for(j in 1:ni){
      p2 <- p2 + deltai[j]*(t(beta)%*%Hi[,j]-log(sum((yi>=yi[j])*t(exp(t(beta)%*%Hi)))))
    }
    p <- p + theta[i]*(p1-alpha*p2)
  }
  p3 <- sum(abs(beta))
  
  loss <- p + p3
  
  return(loss)
}