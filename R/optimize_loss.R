
#' @export
optimize_loss <- function(X,H0=NULL,k,y,delta,theta,alpha,lambda,eta=1,tol=0.01,maxit=5000,tol_H=1e-4,maxit_H=10000,step=1e-5,mu=.99){
  #initialize
  if(is.null(H0)){
    H0 <- init_H(X,k,y,delta,theta,alpha,lambda,eta,method='random',ninit=5)
  }
  H <- H0
  loss <- 0
  eps <- 1
  it <- 0
  
  while(eps > tol & it < maxit){
    loss_prev <- loss
    eps_prev <- eps
    
    # Update beta
    beta <- update_beta(H,y,delta,theta,lambda,eta)
    
    # Update W
    pr <- tryCatch({
      W <- t(update_W(X,H,theta))},
      error=function(e) e)
    if(inherits(pr,"error")){
      break
    }
    
    # Update H
    H <- grad_desc_H(X,W,H,beta,alpha,y,delta,theta,tol_H,maxit_H,step,mu)
    
    # Calculate loss
    l <- calc_loss(X,W,H,beta,alpha,y,delta,theta,lambda)
    loss <- l$loss
    
    eps <- abs(loss - loss_prev)
    
    it <- it + 1
    if(it==maxit){
      warning("Iteration limit reached without convergence")
    }
    
    
    print(sprintf("iter: %d eps: %.4f eps_prev: %.4f",it,eps,eps_prev))
    
  }
  
  return(list(beta=beta,H=H,H0=H0,W=W,loss=loss,eps=eps,survloss=l$surv_loss,recon_err=l$nmf_loss,pen=l$pen_loss))
}