#' @export
update_beta <- function(H,y,delta,theta,lambda){
  ns <- numeric()
  M <- length(H)
  weights <- numeric()
  for(i in 1:M){
    ns[i] <- ncol(H[[i]])
    weights <- c(weights,rep(theta[i],ns[i]))
  }
  H <- t(do.call('cbind',H))
  y <- unlist(y)
  delta <- unlist(delta)
  
  y_surv <- Surv(y,delta)
  fit <- glmnet(H,y_surv,family='cox',weights=weights,lambda=lambda)
  
  return(as.numeric(coef(fit)))
}