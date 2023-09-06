# select highly expressed and variable genes
dat <- matrix(1,nrow=5000,ncol=141)

# consensus nmf
nrun <- 200
k <- 6
ccmat <- matrix(0,nrow=nrow(dat),ncol=nrow(dat))
set.seed(5)

for(rep in 1:nrun){
  cols <- sample(1:ncol(dat),round(.8*ncol(dat)),replace = FALSE)
  keep <- dat[,cols]
  fit1 <- nmf(keep,k,method='lee',nrun=20,maxIter=10,.options=list(verbose=TRUE))
  fit2 <- nmf(keep,k,method='lee',nrun=1,seed=fit1,.options=list(verbose=TRUE))
  W <- fit2@fit@W
  H <- fit2@fit@H
  #pick out exemplar genes from each factor
  sparsity <- matrix(0,nrow=nrow(dat),ncol=k)
  for(j in 1:k){
    sparsity[,j] <- W[,j]-max(W[,setdiff(1:k,j)])
  }
  for(j in 1:k){
    sfilter <- which(sparsity[,j] > quantile(sparsity[,j],1-50/nrow(dat)))
    ccmat[sfilter,sfilter] <- ccmat[sfilter,sfilter] + 1
  }
  print(sprintf('iteration: %d',rep))
}
