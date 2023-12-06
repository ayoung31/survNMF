# source('R/optimize_loss.R')
# source('R/init_H.R')
library(dplyr)
library(NMF)

# let alpha, lambda, k be a grid of possible values

#' @export
cv <- function(X,y,delta,theta,nfold,alpha,lambda=NULL,eta,seed,folds,k){
  library(cvwrapr)
  set.seed(seed)
  M <- length(X)
  loss <- data.frame(matrix(ncol=10,nrow=0))
  colnames(loss) <- c('fold','k','alpha','lambda','loss','nmf_loss','surv_loss','pen_loss','nbeta','cindex')
  r <- 1
  for(f in 1:nfold){
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
    
    for(a in alpha){
      for(l in lambda){
        for(e in eta){
          #fit loss function to training data 
          fit <- optimize_loss(X=Xtrain,k=k,y=ytrain,delta=dtrain,theta=theta,alpha=a,lambda=l,eta=e)
          
          #get Htest
          Htest <- list()
          for(m in 1:M){
            Htest[[m]] <- .fcnnls(fit$W,Xtest[[m]])$coef
          }
          
          #calculate loss for Xtest
          testloss <- calc_loss(Xtest,fit$W,Htest,fit$beta,a,ytest,dtest,theta,l)
          
          #calculate c-index
          Htest_mat <- t(do.call('cbind',Htest))
          ci <- cvwrapr::getCindex(Htest_mat %*% fit$beta, Surv(unlist(ytest),unlist(dtest)))
          
          
          loss[r,] <- c(f,k,a,l,testloss$loss,testloss$nmf_loss,testloss$surv_loss,testloss$pen_loss,sum(fit$beta > 0),ci)
          r <- r+1
          print(sprintf('k: %d fold: %d alpha: %.2f lambda: %.2f eta: %.2f',k,f,a,l,e))
        }
        
      }
    }
  }
  #aloss <- loss %>% group_by(k,alpha,lambda) %>% summarise(avgloss=mean(loss))
  # params <- aloss[which.min(aloss$avgloss),]
  
  return(loss)
}