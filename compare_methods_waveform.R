rm(list=ls())
set.seed(1234)
library(randomForest)
library(sm) ## kde 
source("funcs.R")
source("funcs_patrasen.R")
source("funcs_scott.R")
library(xtable)

################################################ varying feature #####

load('data/waveform.RData')

dat1 = dat1[1:2000,]
dat2 = dat2[1:2000,]

alpha = 0.6

m = floor(min(nrow(dat1),nrow(dat2))/(2-alpha)) ## size of positive samples
n1 = floor(m*(1-alpha))
n2 = m - n1
cl = c(rep(1,n1),rep(0,n2)) ## true class label of unlabeled samples

x0 = rbind(dat1[1:n1,],dat2[1:n2,]) ## unlabeled samples
x1 = dat1[(n1+1):(n1+m),] ## positive samples

dat = rbind(x1,x0)
dat$cl <- as.factor(c(rep("labeled",nrow(x1)), rep("unlabeled",nrow(x0))))
#pihat <- mean(dat$cl=="labeled")

## remove good features to make max(patra/sen) bad
dat = dat[,-c(8,9,15,16,17,18,19)]

nfeat = ncol(dat) - 1
alphas_all = rep(0,nfeat)
ind = (dat$cl == 'labeled')
for(ii in 1:nfeat){
  alphas_all[ii] = method_c_patrasen(dat[!ind,ii],dat[ind,ii],cl,alpha,F)[[1]]
}
plot(alphas_all,ylim=c(0,1))
abline(h=alpha,col=2)
max(alphas_all)

set.seed(2017)
tmpcl = c(rep(1,m),rep(0,n1+n2))
rf.fit <- randomForest(as.factor(tmpcl)~.,data=dat[,1:nfeat], ntree=1000)
imp = importance(rf.fit)
rf.preds <- predict(rf.fit,type='prob')
preds <- rf.preds[,2]
lu <- rep(0,length(preds))
lu[dat$cl=="labeled"] <- 1
p0 = preds[lu==0]
p1 = preds[lu==1]
rst2 = method_c_storey(p0, p1, cl, alpha,T)
rst3 = method_c_patrasen(p0, p1, cl, alpha,T)
rst4 = method_scott(x0[,1:nfeat],x1[,1:nfeat],cl,alpha)
rst5 = method_spy(x0[,1:nfeat], x1[,1:nfeat], cl) ## spy
alpha2 = rst2$alpha
alpha3 = rst3$alpha
alpha4 = rst4$alpha
alpha5 = rst5$alpha

ii = 8
tmp0_8 = dat[!ind,ii]
tmp1_8 = dat[ind,ii]
ii = 5
tmp0_5 = dat[!ind,ii]
tmp1_5 = dat[ind,ii]

save(alphas_all, nfeat, alpha, alpha2, alpha3, alpha4, alpha5, imp, 
     tmp0_8, tmp1_8, tmp0_5,tmp1_5, file = 'data/waveform_vary_feature.RData')

################################################ varying alpha #####

load('data/waveform.RData')

dat1 = dat1[1:2000,]
dat2 = dat2[1:2000,]

alphas = seq(0.01,.99,length.out = 99)

##
library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer

### c-storey + c-patra/sen + spy #####

registerDoParallel(cl)
pt = proc.time()
alpha_est <- foreach(ii = 1:99, .combine = rbind, 
                     .packages = c("randomForest","Iso", "pdist", "LiblineaR","MASS","sm")) %dopar% 
{
  #ii = 40
  alpha = alphas[ii]
  
  m = floor(min(nrow(dat1),nrow(dat2))/(2-alpha))
  n1 = floor(m*(1-alpha))
  n2 = m - n1
  
  x0 = rbind(dat1[1:n1,],dat2[1:n2,])
  x1 = dat1[(n1+1):(n1+m),]
  cl = c(rep(1,n1),rep(0,n2))
  
  dat = rbind(x1,x0)
  dat$cl <- as.factor(c(rep("labeled",nrow(x1)), rep("unlabeled",nrow(x0))))
  pihat <- mean(dat$cl=="labeled")
  
  set.seed(17)
  rf.fit <- randomForest(cl~.,data=dat, ntree=1000)
  rf.preds <- predict(rf.fit,type='prob')
  preds <- rf.preds[,1]
  lu <- rep(0,length(preds))
  lu[dat$cl=="labeled"] <- 1
  p0 = preds[lu==0]
  p1 = preds[lu==1]
  rst3 = method_c_patrasen(p0, p1, cl, alpha, plotit = F) ## patra/sen
  rst5 = method_spy(x0[,-22], x1[,-22], cl) ## spy
  
  rst2 = method_c_storey(p0, p1, cl, alpha, plotit = F) ## storey
 
  #rst4 = method_scott(x0[,-22], x1[,-22], cl, alpha) ## scott
  
  c(alpha, unlist(rst2),unlist(rst3),unlist(rst5))
}
proc.time() - pt ## 416
stopImplicitCluster()

tbs = alpha_est[,4:7]
f1score2s = 2*tbs[,1]/(2*tbs[,1]+tbs[,2]+tbs[,3])
tbs = alpha_est[,10:13]
f1score3s = 2*tbs[,1]/(2*tbs[,1]+tbs[,2]+tbs[,3])
tbs = alpha_est[,16:19]
f1score5s = 2*tbs[,4]/(2*tbs[,4]+tbs[,2]+tbs[,3])

alpha_est_storey_patrasen_spy = cbind(alpha_est[,c(1:3)],f1score2s,
                                      alpha_est[,c(8:9)],f1score3s,
                                      alpha_est[,c(14:15)],f1score5s)
alpha_est_storey_patrasen_spy = as.data.frame(alpha_est_storey_patrasen_spy)
names(alpha_est_storey_patrasen_spy) = c('alpha', 'alpha2', 'acc2', 'f1score2', 'alpha3', 'acc3','f1score3',
                                         'alpha5', 'acc5','f1score5')
save(alpha_est_storey_patrasen_spy, file = 'data/waveform_vary_alpha_storey_patrasen_spy.RData')

### scott #####

registerDoParallel(cl)
pt = proc.time()
alpha_est_scott <- foreach(ii = 1:99, .combine = rbind, 
                     .packages = c("randomForest","Iso", "pdist", "LiblineaR","MASS","sm")) %dopar% 
  {
    alpha = alphas[ii]
    
    m = floor(min(nrow(dat1),nrow(dat2))/(2-alpha))
    n1 = floor(m*(1-alpha))
    n2 = m - n1
    
    kk = 22
    x0 = rbind(dat1[1:n1,-kk],dat2[1:n2,-kk])
    x1 = dat1[(n1+1):(n1+m),-kk]
    cl = c(rep(1,n1),rep(0,n2))
    
    rst4 = method_scott(x0, x1, cl, alpha, nrep = 5) ## scott
    c(alpha, unlist(rst4))
}
proc.time() - pt ## 
stopImplicitCluster()

tbs = alpha_est_scott[,4:7]
f1score4s = 2*tbs[,1]/(2*tbs[,1]+tbs[,2]+tbs[,3])

alpha_est_scott = cbind(alpha_est_scott[,c(1:3)],f1score4s,alpha_est_scott[,c(4:7)])
alpha_est_scott = as.data.frame(alpha_est_scott)
names(alpha_est_scott) = c('alpha', 'alpha4', 'acc4', 'f1score4','v1','v2','v3','v4')
save(alpha_est_scott, file = 'data/waveform_vary_alpha_scott.RData')

################################################ varying sample size #####

load('data/waveform.RData')

## 
alphas = c(.1,.5,.9)
ns = ceiling(50*2^c(1:7))

#alphas = c(.9)
#ns = ceiling(50*2^c(7))

##
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
cl <- makeCluster(cores[1]-1) #not to overload your computer

## C-storey and C-patra/sen #####

registerDoParallel(cl)
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
  cl = c(rep(1,m1),rep(0,m2))
  
  dat = rbind(x1,x0)
  dat$cl <- as.factor(c(rep("labeled",nrow(x1)), rep("unlabeled",nrow(x0))))
  pihat <- mean(dat$cl=="labeled")
  
  set.seed(17)
  rf.fit <- randomForest(cl~.,data=dat)
  rf.preds <- predict(rf.fit,type='prob')
  preds <- rf.preds[,1]
  lu <- rep(0,length(preds))
  lu[dat$cl=="labeled"] <- 1
  p0 = preds[lu==0]
  p1 = preds[lu==1]
  
  rst2 = method_c_storey(p0, p1, cl, alpha, plotit = F) ## storey
  rst3 = method_c_patrasen(p0, p1, cl, alpha, plotit = F) ## patra/sen
  #rst4 = method_scott(x0[,-22], x1[,-22], cl, alpha) ## scott
  #rst5 = method_spy(x0[,-22], x1[,-22], cl) ## spy
  c(alpha, n0, grp, unlist(rst2),unlist(rst3))
}
proc.time() - pt ## 737
stopImplicitCluster()

alpha_est_patra_storey = as.data.frame(alpha_est)
names(alpha_est_patra_storey) = c('alpha', 'n0', 'grp', 
                                  'alpha2','acc2','v1','v2','v3','v4',
                                  'alpha3','acc3','v1','v2','v3','v4')
save(alpha_est_patra_storey, file = 'data/waveform_vary_size_patra_storey.RData')

## spy #####

registerDoParallel(cl)
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
  cl = c(rep(1,m1),rep(0,m2))
  
  dat = rbind(x1,x0)
  dat$cl <- as.factor(c(rep("labeled",nrow(x1)), rep("unlabeled",nrow(x0))))
  pihat <- mean(dat$cl=="labeled")
  
  set.seed(17)
  rf.fit <- randomForest(cl~.,data=dat)
  rf.preds <- predict(rf.fit,type='prob')
  preds <- rf.preds[,1]
  lu <- rep(0,length(preds))
  lu[dat$cl=="labeled"] <- 1
  p0 = preds[lu==0]
  p1 = preds[lu==1]
  
  #rst2 = method_c_storey(p0, p1, cl, alpha, plotit = F) ## storey
  #rst3 = method_c_patrasen(p0, p1, cl, alpha, plotit = F) ## patra/sen
  #rst4 = method_scott(x0[,-22], x1[,-22], cl, alpha) ## scott
  rst5 = method_spy(x0[,-22], x1[,-22], cl) ## spy
  c(alpha, n0, grp, unlist(rst5))
}
proc.time() - pt ## 1361.755
stopImplicitCluster()

alpha_est_spy = as.data.frame(alpha_est_spy)
names(alpha_est_spy) = c('alpha', 'n0', 'grp', 'alpha5', 'acc5','v1','v2','v3','v4')
save(alpha_est_spy, file = 'data/waveform_vary_size_spy.RData')

## scott #####

registerDoParallel(cl)
pt = proc.time()
set.seed(17)
alpha_est_scott <- foreach(ii = 1:length(index), .combine = rbind, .packages = c("MASS","pdist","LiblineaR")) %dopar% {
  alpha = index[[ii]][1]
  n0 = index[[ii]][2]
  n1 = index[[ii]][3]
  m1 = index[[ii]][4]
  m2 = index[[ii]][5]
  grp = index[[ii]][6]

  x1 = dat1[index[[ii]][7:(6+n1)],]
  x0 = rbind(dat1[index[[ii]][(6+n1+1):(6+n1+m1)],], dat2[index[[ii]][((6+n1+m1+1):((6+n1+m1+m2)))],])

  rst4 = method_scott(x0[,-22], x1[,-22], alpha, nrep = 1)
  c(alpha, n0, grp, unlist(rst4))
}
proc.time() - pt ## 1487
stopImplicitCluster()
alpha_est_scott = as.data.frame(alpha_est_scott)
names(alpha_est_scott) = c('alpha', 'n0', 'grp', 'alpha4','acc4','v1','v2','v3','v4')
save(alpha_est_scott, file = 'data/waveform_vary_size_scott.RData')




