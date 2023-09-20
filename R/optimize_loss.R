library(NMF)
library(survival)
library(glmnet)

# source('R/update_beta.R')
# source('R/update_W.R')
# source('R/gradient_descent_H.R')
# source('R/calc_loss.R')
# source('R/init_H.R')

#' @export
optimize_loss <- function(X,H0,k,y,delta,theta,alpha,lambda,tol=0.001,maxit=5000,tol_H=1e-4,maxit_H=15000,step=1e-6,mu=.9){
  #initialize
  H <- H0
  loss <- 0
  eps <- 1
  it <- 0
  
  while(eps > tol & it < maxit){
    loss_prev <- loss
    
    # Update beta
    beta <- update_beta(H,y,delta,theta,lambda)
    
    # Update W
    W <- t(update_W(X,H,theta))
    
    # Update H
    H <- grad_desc_H(X,W,H,beta,alpha,y,delta,theta,tol_H,maxit_H,step,mu)
    
    # Calculate loss
    loss <- calc_loss(X,W,H,beta,alpha,y,delta,theta,lambda)
    
    eps <- abs(loss - loss_prev)
    
    it <- it + 1
    if(it==maxit){
      warning("Iteration limit reached without convergence")
    }
    print(sprintf("iter: %d eps: %.3f",it,eps))
  }
  
  return(list(beta=beta,H=H,W=W,loss=loss,eps=eps))
}