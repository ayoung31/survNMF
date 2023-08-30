source('R/optimize_loss.R')
source('R/sim_dat.R')

set.seed(35)

M <- 3
p <- 1000
ni <- 100
k <- 6
ksurv <- 2
beta <- rep(1,ksurv)

dat <- sim_dat(M,p,ni,k,ksurv,beta)

X <- dat$X
y <- dat$y
delta <- dat$delta

start <- Sys.time()
fit <- optimize_loss(X,k,alpha,y,delta,theta,lambda,tol,maxit,tol_H,maxit_H,step)
time <- Sys.time() - start

save(dat,fit,time,file='data/test_optim.RData')