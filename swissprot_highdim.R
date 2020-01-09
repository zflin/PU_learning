########################################################################
# tcdb_swissprot
########################################################################

## goal: compare C-PS and PS on swissprot data with
## different number of simulated noise features and
## different amounts of variable screening

rm(list=ls())
set.seed(1234)
library(randomForest)
source("funcs.R")
source("funcs_roc.R")
source("funcs_patrasen.R")

load("data/tcdb_swissprot.RData")

y <- dat$cl == 'labeled'
feats <- as.matrix(dat[,!(colnames(dat)=="cl")])
kvar <- c(500,200,100,50,10,1)
n_append <- c(-1,ncol(dat)-1,9*(ncol(dat)-1))
res <- matrix(0,ncol=2*length(n_append),nrow=length(kvar))

AppendNoisyFeatures <- function(feats,nfake){
  noisy_features <- matrix(NA_real_,nrow=nrow(feats),ncol=nfake)
  col_samples <- sample(1:ncol(feats),nfake,replace=TRUE)
  for(ii in 1:ncol(noisy_features)){
    noisy_features[,ii] <- sample(feats[,col_samples[ii]])
  }
  return(cbind(feats,noisy_features))
}

p_vals_all <- list()
for(ii in 1:length(n_append)){
  print(paste0("n_append:",n_append[ii]))
  ## append features
  if(n_append[ii]>0){
    feats <- AppendNoisyFeatures(feats,n_append[ii])
  }
  ## compute p-values for each feature
  ## use fisher exact or chisq
  pvals <- rep(-1,ncol(feats))
  for(kk in 1:length(pvals)){
    if(min(table(feats[,kk],y)) < 10){
      pvals[kk] <- fisher.test(feats[,kk],y)$p.value  
    } else {
      pvals[kk] <- chisq.test(feats[,kk],y)$p.value
    }
  }
  p_vals_all[[ii]] <- pvals
  ## run for different amounts of tuning
  for(jj in 1:nrow(res)){
    print(paste0("n_features:",kvar[jj]))
    ## C-PS
    feats_temp <- feats[,order(pvals)[1:kvar[jj]],drop=FALSE]
    rffit <- randomForest(feats_temp,as.factor(y))
    preds <- predict(rffit,type="prob")[,1]
    rst3 <- method_c_patrasen(preds[!y], preds[y], cl)
    ## max single feature PS 
    alphas_all = rep(NA_real_,ncol(feats_temp))
    for(kk in 1:ncol(feats_temp)){
      alphas_all[kk] <- method_c_patrasen(feats_temp[!y,kk],feats_temp[y,kk])[[1]]
    }
    ## store results
    res[jj,2*(ii-1)+1] <- rst3$alpha
    res[jj,2*(ii-1)+2] <- max(alphas_all)
  }
}
rownames(res) <- kvar
colnames(res) <- rep(c("C-PS","PS"),3)
save(res,n_append,kvar,p_vals_all,file="results/swissprot_highdim.RData")