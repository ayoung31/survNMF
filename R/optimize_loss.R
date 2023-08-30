library(NMF)
library(survival)
library(glmnet)

source('R/update_beta.R')
source('R/update_W.R')
source('R/gradient_descent_H.R')
source('R/calc_loss.R')
source('R/init_H.R')

optimize_loss <- function(X,k,alpha,y,delta,theta,lambda,tol,maxit,tol_H,maxit_H,step){
  #initialize
  H <- init_H(X,k)
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
    H <- grad_desc_H(X,W,H,beta,alpha,y,delta,theta,tol_H,maxit_H,step)
    
    # Calculate loss
    loss <- calc_loss(X,W,H,beta,alpha,y,delta,theta,lambda)
    
    eps <- abs(loss - loss_prev)
    
    it <- it + 1
    if(it==maxit){
      warning("Iteration limit reached without convergence")
    }
    print(sprintf("iter: %d eps: %.3f",it,eps))
  }
  
  return(list(beta,H,W,loss,eps))
}