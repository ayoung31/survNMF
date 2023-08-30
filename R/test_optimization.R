source('R/optimize_loss.R')
source('R/sim_dat.R')

set.seed(35)

M <- 3
p <- 500
ni <- 100
k <- 6
ksurv <- 2
beta <- rep(1,ksurv)

dat <- sim_dat(M,p,ni,k,ksurv,beta)

X <- dat$X
y <- dat$y
delta <- dat$delta

start <- Sys.time()
print(start)
fit <- optimize_loss(X,k,alpha=1,y,delta,theta=c(1,1,1),lambda=0,tol=0.01,maxit=5000,tol_H=1e-4,maxit_H=15000,step=1e-6)
time <- Sys.time() - start

save(dat,fit,time,file='data/test_optim.RData')