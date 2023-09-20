source('R/cross_validation.R')

load('data/sim_dat.RData')

test <- cv(X,y,delta,theta=c(1,1,1),nfold=5,alpha=1,K=4:8,seed=123)

save(test,file='data/cv_res.RData')