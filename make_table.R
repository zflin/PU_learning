### produce output table for real data
### and summaries of results
rm(list=ls())
library(xtable)
source("funcs.R")
source("funcs_patrasen.R")
load("results/swissprot_highdim.RData")
load("data/tcdb_swissprot.RData")

#### MAKE TABLE
colnames(res)[colnames(res)=="max(PS)"] <- "PS"

res <- data.frame(k=rownames(res),res,check.names=FALSE)
out <- xtable(res,
              caption="$\\alpha$ estimates from C-PS and single feature PS with p, 2p, and 10p features when prescreening the top k features. True $\\alpha \\approx 0.93$.",
              align="cc|cc|cc|cc",
              label="tab:tcdb_swissprot")
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- c("& \\multicolumn{2}{c}{p} & \\multicolumn{2}{c}{2p} & \\multicolumn{2}{c}{10p} \\\\\n")
print(out,
      file="figs/comparison_varying_data.tex",
      add.to.row=addtorow,
      include.rownames=FALSE)

#### FEATURE COMMONALITIES
for(k in kvar){
  num_common <- sum(order(p_vals_all[[1]])[1:k] %in% order(p_vals_all[[3]])[1:k])
  print(paste0(num_common,"/",k," common features"))
}


## compute 2 x 2 tables, p-values, PS alphas for
## a few features
y <- dat$cl == 'labeled'
feats <- as.matrix(dat[,!(colnames(dat)=="cl")])
## print 2 x 2 table for features
## lowest p-value and highest alpha_hat

alpha_hats <- rep(NA_real_,ncol(feats))
for(ii in 1:ncol(feats)){
  alpha_hats[ii] <- method_c_patrasen(feats[!y,ii],feats[y,ii])[[1]]
}

## smallest p feature
ix <- which.min(p_vals_all[[1]])
table(feats[,ix],y)
alpha_hats[ix]
p_vals_all[[1]][ix]

## largest alpha feature in top 50
ix <- which.max(1*(rank(p_vals_all[[1]])<=10)*alpha_hats)
table(feats[,ix],y)
alpha_hats[ix]
p_vals_all[[1]][ix]
sum(p_vals_all[[1]] < p_vals_all[[1]][ix])


## check on method_c_patrasen for binary features
BinaryPS <- function(p0,p1){
  return(1-min(mean(p0==0)/mean(p1==0),mean(p0==1)/mean(p1==1)))
}
# alpha_hats_binary <- rep(NA_real_,ncol(feats))
# for(ii in 1:ncol(feats)){
#   alpha_hats_binary[ii] <- BinaryPS(feats[!y,ii],feats[y,ii])[[1]]
# }
# plot(alpha_hats_binary,alpha_hats)
