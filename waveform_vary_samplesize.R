########################################################################
# varying sample size
# Note: This is a time-consuming Rscript, which takes about 6000 seconds
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

## alpha and sample size setting
alphas = c(.1,.5,.9)
ns = ceiling(50*2^c(1:7))

## index contains setting info of each round of calculation
index = list()
kk = 1
grp = 1
set.seed(17)
for(ii in 1:length(alphas)){
  for(jj in 1:length(ns)){
    alpha = alphas[ii] ## proportion of 2 in unlabeled set
    n0 = ns[jj] ## size of unlabeled set
    m1 = ceiling(n0*(1-alpha)) ## number of 1 in unlabeled set
    m2 = n0 - m1 ## number of 2 in unlabeled set
    n1 = n0 ## size of labeled set
    for(ll in 1:20){
      index[[kk]] = c(alpha, n0, n1, m1, m2, grp, sample(nrow(dat1), n1+m1), sample(nrow(dat2), m2))
      kk = kk + 1
    }
    grp = grp + 1
  }
}

##
library(foreach)
library(doParallel)
cores=detectCores()
cls <- makeCluster(cores[1]-1) #not to overload your computer

## C-roc and C-patra/sen #####
set.seed(17)

registerDoParallel(cls)
pt = proc.time()
alpha_est <- foreach(ii = 1:length(index), .combine = rbind, .packages = c("randomForest","Iso")) %dopar% {
  alpha = index[[ii]][1]
  n0 = index[[ii]][2]
  n1 = index[[ii]][3]
  m1 = index[[ii]][4]
  m2 = index[[ii]][5]
  grp = index[[ii]][6]
  
  x1 = dat1[index[[ii]][7:(6+n1)],]
  x0 = rbind(dat1[index[[ii]][(6+n1+1):(6+n1+m1)],], dat2[index[[ii]][((6+n1+m1+1):((6+n1+m1+m2)))],])
  
  dat = rbind(x1,x0)
  dat$cl <- as.factor(c(rep("labeled",nrow(x1)), rep("unlabeled",nrow(x0))))
  pihat <- mean(dat$cl=="labeled")
  
  rf.fit <- randomForest(cl~.,data=dat)
  rf.preds <- predict(rf.fit,type='prob')
  preds <- rf.preds[,1]
  lu <- rep(0,length(preds))
  lu[dat$cl=="labeled"] <- 1
  p0 = preds[lu==0]
  p1 = preds[lu==1]
  
  rst2 = method_c_roc(p0, p1) ## storey
  rst3 = method_c_patrasen(p0, p1) ## patra/sen
  c(alpha, n0, grp, unlist(rst2),unlist(rst3))
}
proc.time() - pt ## 537.217
stopImplicitCluster()

alpha_est_croc_patrasen = as.data.frame(alpha_est)
names(alpha_est_croc_patrasen) = c('alpha', 'n0', 'grp', 
                                  'alpha2','acc2','tb2',
                                  'alpha3','acc3','tb3')
save(alpha_est_croc_patrasen, file = 'results/waveform_vary_size_croc_patrasen.RData')

## spy #####
set.seed(17)

registerDoParallel(cls)
pt = proc.time()
alpha_est_spy <- foreach(ii = 1:length(index), .combine = rbind, .packages = c("randomForest")) %dopar% {
  alpha = index[[ii]][1]
  n0 = index[[ii]][2]
  n1 = index[[ii]][3]
  m1 = index[[ii]][4]
  m2 = index[[ii]][5]
  grp = index[[ii]][6]

  x1 = dat1[index[[ii]][7:(6+n1)],]
  x0 = rbind(dat1[index[[ii]][(6+n1+1):(6+n1+m1)],], dat2[index[[ii]][((6+n1+m1+1):((6+n1+m1+m2)))],])

  rst5 = method_spy(x0[,-22], x1[,-22]) ## spy
  c(alpha, n0, grp, unlist(rst5))
}
proc.time() - pt ## 732.968
stopImplicitCluster()

alpha_est_spy = as.data.frame(alpha_est_spy)
names(alpha_est_spy) = c('alpha', 'n0', 'grp', 'alpha5', 'acc5','tb5')
save(alpha_est_spy, file = 'results/waveform_vary_size_spy.RData')

## roc #####
set.seed(17)

registerDoParallel(cls)
pt = proc.time()
alpha_est_roc <- foreach(ii = 1:length(index), .combine = rbind, .packages = c("MASS","pdist","LiblineaR")) %dopar% {
  alpha = index[[ii]][1]
  n0 = index[[ii]][2]
  n1 = index[[ii]][3]
  m1 = index[[ii]][4]
  m2 = index[[ii]][5]
  grp = index[[ii]][6]

  x1 = dat1[index[[ii]][7:(6+n1)],]
  x0 = rbind(dat1[index[[ii]][(6+n1+1):(6+n1+m1)],], dat2[index[[ii]][((6+n1+m1+1):((6+n1+m1+m2)))],])

  rst4 = method_roc(x0[,-22], x1[,-22])
  c(alpha, n0, grp, unlist(rst4))
}
proc.time() - pt ## 3956.811 seconds
stopImplicitCluster()
alpha_est_roc = as.data.frame(alpha_est_roc)
names(alpha_est_roc) = c('alpha', 'n0', 'grp', 'alpha4','acc4','tb4')
save(alpha_est_roc, file = 'results/waveform_vary_size_roc.RData')




