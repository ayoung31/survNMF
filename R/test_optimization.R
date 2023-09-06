source('R/optimize_loss.R')
source('R/sim_dat.R')
source('R/init_H.R')

set.seed(35)

M <- 3
p <- 500
ni <- 100
k <- 6
ksurv <- 2
beta <- rep(2,ksurv)

dat <- sim_dat(M,p,ni,k,ksurv,beta)

X <- dat$X
y <- dat$y
delta <- dat$delta

H0 <- init_H(X,k)

for(lambda in seq(0,.4,.02)){
  start <- Sys.time()
  print(start)
  fit <- optimize_loss(X,H0,k,y,delta,theta=c(1,1,1),alpha=1,lambda=lambda,tol=0.001,maxit=5000,tol_H=1e-4,maxit_H=15000,step=1e-6)
  time <- Sys.time() - start
  print(time)
  
  save(dat,fit,time,file=paste0('data/test_optim',lambda,'.RData'))
}


