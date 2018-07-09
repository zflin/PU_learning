#######################################################################
# varying feature
########################################################################

rm(list=ls())
set.seed(1234)
library(randomForest)
source("funcs.R")
source("funcs_patrasen.R")
source("funcs_roc.R")
library(xtable)

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

## remove good features to make max(patra/sen) bad
dat = dat[,-c(8,9,15,16,17,18,19)]

## use patra/sen on each feature
nfeat = ncol(dat) - 1
alphas_all = rep(0,nfeat)
ind = (dat$cl == 'labeled')
for(ii in 1:nfeat){
  alphas_all[ii] = method_c_patrasen(dat[!ind,ii],dat[ind,ii],cl,alpha,F)[[1]]
}
# plot(alphas_all,ylim=c(0,1))
# abline(h=alpha,col=2)
# max(alphas_all) ## naive patra/sen estimator

## use C-patra/sen, C-roc, SPY, ROC on entire data
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
rst2 = method_c_roc(p0, p1, cl, alpha)
rst3 = method_c_patrasen(p0, p1, cl, alpha)
rst4 = method_roc(x0[,1:nfeat],x1[,1:nfeat],cl,alpha)
rst5 = method_spy(x0[,1:nfeat], x1[,1:nfeat], cl) ## spy
alpha2 = rst2$alpha
alpha3 = rst3$alpha
alpha4 = rst4$alpha
alpha5 = rst5$alpha

## save 8th and 5th features
ii = 8
tmp0_8 = dat[!ind,ii]
tmp1_8 = dat[ind,ii]
ii = 5
tmp0_5 = dat[!ind,ii]
tmp1_5 = dat[ind,ii]

save(alphas_all, nfeat, alpha, alpha2, alpha3, alpha4, alpha5, imp, 
     tmp0_8, tmp1_8, tmp0_5,tmp1_5, file = 'data/waveform_vary_feature.RData')
