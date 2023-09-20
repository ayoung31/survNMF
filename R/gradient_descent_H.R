# This function performs gradient descent on one column of H

#' @export
f_grad <- function(X,W,H,beta,alpha,y,delta,theta,j){
  # X is the p x n_i expression matrix of study i
  # W is a matrix of p genes by k factors
  # H is a matrix k factors by n_i subjects for study i
  # beta is a k x 1 vector of regression coefficient
  # alpha is a tuning parameter representing contribution of survival to loss fctn
  # y is a n_i x 1 vector of censored survival outcomes
  # delta is a n_i x 1 vector of indicators, 1=event, 0=censored
  # theta is the specified weight for study i
  # j is the current subject
  
  p1 <- t(W)%*%W%*%H[,j] - t(W)%*%X[,j]
  p2 <- as.numeric(1 - (exp(t(beta)%*%H[,j]))/(sum((y>=y[j])*t(exp(t(beta)%*%H)))))
  
  return(theta*(p1-alpha*delta[j]*p2*beta))
}

#' @export
grad_desc_hji <- function(X,W,H,beta,alpha,y,delta,theta,j,tol,maxit,step,mu){
  it <- 0
  start <- Sys.time()
  eps <- 1
  b <- 0
  while(eps > tol & it < maxit){
    grad <- f_grad_rcpp(X,W,H,as.matrix(beta),alpha,as.matrix(y),as.matrix(delta),theta,j)
    #grad <- f_grad(X,W,H,beta,alpha,y,delta,theta,j)
    Hj_prior <- H[,j]
    b_prior <- b
    b <- mu*b_prior + grad
    H[,j] <- H[,j]-step*b
    H[,j][H[,j]<0] <- 0
    
    
    eps <- sqrt(sum((Hj_prior-H[,j])^2))
    if(is.na(eps)){
      print(sprintf('it: %d',it))
    }
    
    it <- it + 1
    if(it==maxit){
      warning("Iteration limit reached without convergence")
    }
    
    if(it %% 50 == 0){
      sprintf("iter: %d eps: %.3f",it,eps)
    }
  }
  run_time <- Sys.time() - start
  
  return(list(H=H,run_time=run_time))
}

# here H,X,y,and delta are lists
# theta is a vector

#' @export
grad_desc_H <- function(X,W,H,beta,alpha,y,delta,theta,tol,maxit,step,mu){
  M <- length(X)
  for(i in 1:M){
    ni <- ncol(X[[i]])
    for(j in 1:ni){
      fit <- grad_desc_hji(X[[i]],W,H[[i]],beta,alpha,y[[i]],delta[[i]],theta[i],j,tol,maxit,step,mu)
      H[[i]] <- fit$H
      sprintf("i: %d j: %d runtime: %d", i, j, fit$runtime)
    }
  }
  return(H)
}