########################################################################
# varying alpha
# Note: This is a time-consuming Rscript, which takes about 2000 seconds
#       with 7 processors on an i7 CPU computer.
#       Please reserve about 10 GB memory for this code.
########################################################################

rm(list=ls())
set.seed(1234)
library(randomForest)
source("funcs.R")
source("funcs_patrasen.R")
source("funcs_roc.R")
library(xtable)

load('data/waveform.RData')

alphas = seq(0.01,.99,length.out = 99)

library(foreach)
library(doParallel)
cores=detectCores()
cls <- makeCluster(cores[1]-1) #not to overload your computer

### c-roc + c-patra/sen #####
set.seed(17)

registerDoParallel(cls)
pt = proc.time()
alpha_est <- foreach(ii = 1:99, .combine = rbind, .packages = c("randomForest","Iso")) %dopar% {

  alpha = alphas[ii]

  n = 3000 ## sample size for U
  m = 3000 ## sample size for P

  n1 = floor(n*(1-alpha))
  n2 = n - n1

  x0 = rbind(dat1[1:n1,],dat2[1:n2,])
  x1 = dat1[(n1+1):(n1+m),]

  dat = rbind(x1,x0)
  dat$cl <- as.factor(c(rep("labeled",nrow(x1)), rep("unlabeled",nrow(x0))))
  pihat <- mean(dat$cl=="labeled")

  rf.fit <- randomForest(cl~.,data=dat, ntree=1000)
  rf.preds <- predict(rf.fit,type='prob')
  preds <- rf.preds[,1]
  lu <- rep(0,length(preds))
  lu[dat$cl=="labeled"] <- 1
  p0 = preds[lu==0]
  p1 = preds[lu==1]

  cl = c(rep(1,n1),rep(0,n2)) ## true labels in U

  rst3 = method_c_patrasen(p0, p1, cl) ## c-patra/sen
  rst2 = method_c_roc(p0, p1, cl) ## c-roc

  c(alpha, unlist(rst2),unlist(rst3))
}
proc.time() - pt ## 449.121 seconds
stopImplicitCluster()

## calculate F1-score
tbs = alpha_est[,4:7]
f1score2s = 2*tbs[,1]/(2*tbs[,1]+tbs[,2]+tbs[,3])
tbs = alpha_est[,10:13]
f1score3s = 2*tbs[,1]/(2*tbs[,1]+tbs[,2]+tbs[,3])

alpha_est_croc_patrasen = cbind(alpha_est[,c(1:3)],f1score2s,
                                      alpha_est[,c(8:9)],f1score3s)
alpha_est_croc_patrasen = as.data.frame(alpha_est_croc_patrasen)
names(alpha_est_croc_patrasen) = c('alpha', 'alpha2', 'acc2', 'f1score2', 'alpha3', 'acc3','f1score3')
save(alpha_est_croc_patrasen, file = 'data/waveform_vary_alpha_croc_patrasen.RData')

### spy #####
set.seed(17)

registerDoParallel(cls)
pt = proc.time()
alpha_est <- foreach(ii = 1:99, .combine = rbind, .packages = c("randomForest")) %dopar% {
  
  alpha = alphas[ii]
  
  n = 3000 ## sample size for U
  m = 3000 ## sample size for P
  
  n1 = floor(n*(1-alpha))
  n2 = n - n1
  
  x0 = rbind(dat1[1:n1,],dat2[1:n2,])
  x1 = dat1[(n1+1):(n1+m),]
  
  cl = c(rep(1,n1),rep(0,n2)) ## true labels in U
  
  rst5 = method_spy(x0[,-22], x1[,-22], cl) ## spy
  
  c(alpha, unlist(rst5))
}
proc.time() - pt ## 231.489 seconds
stopImplicitCluster()

## calculate F1-score
tbs = alpha_est[,4:7]
f1score5s = 2*tbs[,4]/(2*tbs[,4]+tbs[,2]+tbs[,3])

alpha_est_spy = cbind(alpha_est[,c(1:3)], f1score5s)
alpha_est_spy = as.data.frame(alpha_est_spy)
names(alpha_est_spy) = c('alpha', 'alpha5', 'acc5','f1score5')
save(alpha_est_spy, file = 'data/waveform_vary_alpha_spy.RData')

### roc #####
set.seed(17)

registerDoParallel(cls)
pt = proc.time()
alpha_est_roc <- foreach(ii = 1:99, .combine = rbind, .packages = c("pdist", "LiblineaR","MASS")) %dopar% {
  alpha = alphas[ii]

  n = 3000
  m = 3000

  n1 = floor(n*(1-alpha))
  n2 = n - n1

  kk = 22
  x0 = rbind(dat1[1:n1,-kk],dat2[1:n2,-kk])
  x1 = dat1[(n1+1):(n1+m),-kk]
  cl = c(rep(1,n1),rep(0,n2))

  rst4 = method_roc(x0, x1, cl) ## ked
  c(alpha, unlist(rst4))
}
proc.time() - pt ## 2889.328 seconds
stopImplicitCluster()

tbs = alpha_est_roc[,4:7]
f1score4s = 2*tbs[,1]/(2*tbs[,1]+tbs[,2]+tbs[,3])

alpha_est_roc = cbind(alpha_est_roc[,c(1:3)],f1score4s,alpha_est_roc[,c(4:7)])
alpha_est_roc = as.data.frame(alpha_est_roc)
names(alpha_est_roc) = c('alpha', 'alpha4', 'acc4', 'f1score4','v1','v2','v3','v4')
save(alpha_est_roc, file = 'data/waveform_vary_alpha_roc.RData')
