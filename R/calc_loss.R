calc_loss <- function(X,W,H,beta,alpha,y,delta,theta,lambda){
  M <- length(X)
  p1 <- 0
  p2 <- 0
  for(i in 1:M){
    ni <- ncol(X[[i]])
    yi <- y[[i]]
    Xi <- X[[i]]
    Hi <- H[[i]]
    deltai <- delta[[i]]
    p1a <- norm(Xi-W%*%Hi,'F')^2
    p2a <- 0
    for(j in 1:ni){
      p2a <- p2a + deltai[j]*(t(beta)%*%Hi[,j]-log(sum((yi>=yi[j])*t(exp(t(beta)%*%Hi)))))
    }
    p1 <- p1 + theta[i]*p1a
    p2 <- p2 - theta[i]*alpha*p2a
  }
  p3 <- sum(abs(beta))
  
  loss <- p1 + p2 + p3
  survloss <- p2 + p3
  
  return(list(loss=loss,survperc=100*survloss/loss))
}