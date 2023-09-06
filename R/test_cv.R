
#try out cv
source('R/sim_dat.R')
source('R/cross_validation.R')


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
#seq(0,.3,.02)
grid <- cv(X,y,delta,theta=c(1,1,1),nfold=5,alpha=1,lambda=seq(.06,.2,.02),K=4:8,seed=6)

save(grid,file='data/cv_grid.RData')