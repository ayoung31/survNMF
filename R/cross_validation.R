# source('R/optimize_loss.R')
# source('R/init_H.R')
library(dplyr)
library(NMF)

# let alpha, lambda, k be a grid of possible values

#' @export
cv <- function(X,y,delta,theta,nfold,alpha,lambda=NULL,K,seed,folds,f,k){
  M <- length(X)
  loss <- data.frame(matrix(ncol=5,nrow=0))
  colnames(loss) <- c('fold','k','alpha','lambda','loss')
  r <- 1
  Xtrain <- list()
  Xtest <- list()
  ytrain <- list()
  ytest <- list()
  dtrain <- list()
  dtest <- list()
  for(m in 1:M){
    Xtrain[[m]] <- X[[m]][,folds[[m]] != f]
    Xtest[[m]] <- X[[m]][,folds[[m]] == f]
    ytrain[[m]] <- y[[m]][folds[[m]] != f]
    ytest[[m]] <- y[[m]][folds[[m]] == f]
    dtrain[[m]] <- delta[[m]][folds[[m]] != f]
    dtest[[m]] <- delta[[m]][folds[[m]] == f]
  }
  H0 <- init_H(Xtrain,k)
  for(a in alpha){
    for(l in lambda){
      #fit loss function to training data 
      fit <- optimize_loss(X=Xtrain,H0=H0,k=k,y=ytrain,delta=dtrain,theta=theta,alpha=a,lambda=l,tol=0.01,maxit=5000,tol_H=1e-4,maxit_H=10000,step=1e-5,mu=.99)
      
      #get Htest
      Htest <- list()
      for(m in 1:M){
        Htest[[m]] <- .fcnnls(fit$W,Xtest[[m]])$coef
      }
      
      #calculate loss for Xtest
      testloss <- calc_loss(Xtest,fit$W,Htest,fit$beta,a,ytest,dtest,theta,l)
      
      loss[r,] <- c(f,k,a,l,testloss)
      r <- r+1
      print(sprintf('k: %d lambda %.2f',k,l))
    }
  }
  #aloss <- loss %>% group_by(k,alpha,lambda) %>% summarise(avgloss=mean(loss))
  # params <- aloss[which.min(aloss$avgloss),]
  
  return(loss)
}