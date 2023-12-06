#' @export
update_beta <- function(H,y,delta,theta,lambda,eta,stdize=TRUE){
  ns <- numeric()
  M <- length(H)
  weights <- numeric()
  for(i in 1:M){
    ns[i] <- ncol(H[[i]])
    weights <- c(weights,rep(theta[i],ns[i]))
  }
  H <- t(do.call('cbind',H))
  if(!stdize){
    H <- apply(H,2,scale)
  }
  y <- unlist(y)
  delta <- unlist(delta)
  
  y_surv <- Surv(y,delta)
  fit <- glmnet::glmnet(H,y_surv,family='cox',weights=weights,
                        lambda=lambda,alpha=eta,standardize = stdize)
  
  if(!stdize){
    return(list(beta=as.numeric(coef(fit)),H=H))
  }else{
    return(as.numeric(coef(fit)))
  }
}