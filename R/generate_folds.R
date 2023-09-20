#' @export
generate_folds <- function(seed,X){
  set.seed(seed)
  M <- length(X)
  folds <- list()
  for(m in 1:M){
    folds[[m]] <- sample(1:nfold,ncol(X[[m]]),replace=TRUE)
  }
  return(folds)
}

