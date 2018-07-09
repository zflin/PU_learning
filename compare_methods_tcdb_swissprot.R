if(!dir.exists("figs")){
  dir.create("figs")
}

########################################################################
# tcdb_swissprot
########################################################################

rm(list=ls())
set.seed(1234)
library(randomForest)
source("funcs.R")
source("funcs_roc.R")
source("funcs_patrasen.R")
library(xtable)

load("data/tcdb_swissprot.RData")

n1 = 2453
n2 = 4558
n3 = 348
n = n1 + n2 + n3

x1 = dat[1:n1,-ncol(dat)]
x0 = dat[(n1+1):n,-ncol(dat)]

## C-roc, C-patra/sen, ROC, SPY estimators
set.seed(2017)
rf.fit = randomForest(cl~.,data = dat)
rf.preds <- predict(rf.fit,type='prob')
preds <- rf.preds[,1]
lu <- rep(0,length(preds))
lu[dat$cl=="labeled"] <- 1
p0 = preds[lu==0]
p1 = preds[lu==1]

rst2 = method_c_roc(p0, p1, cl, alpha)
rst3 = method_c_patrasen(p0, p1, cl, alpha, plotit = F)
rst4 = method_roc(x0, x1, cl, alpha)
rst5 = method_spy(x0, x1, cl)

## accuracy on all data, treat all positive samples are correctly classified
acc2 = (rst2$acc*(n2+n3) + n1)/n
acc3 = (rst3$acc*(n2+n3) + n1)/n
acc4 = (rst4$acc*(n2+n3) + n1)/n
acc5 = (rst5$acc*(n2+n3) + n1)/n

f1score2 = 2*(rst2$tb[4]+n1)/(rst2$tb[3]+rst2$tb[2]+2*(rst2$tb[4]+n1))
f1score3 = 2*(rst3$tb[4]+n1)/(rst3$tb[3]+rst3$tb[2]+2*(rst3$tb[4]+n1))
f1score4 = 2*(rst4$tb[4]+n1)/(rst4$tb[3]+rst4$tb[2]+2*(rst4$tb[4]+n1))
f1score5 = 2*(rst5$tb[1]+n1)/(rst5$tb[2]+rst5$tb[3]+2*(rst5$tb[1]+n1))

## ideal case: 10 fold cross-validation
dat2 = dat
dat2$cl = as.factor(c(rep("labeled",n1),rep("unlabeled",n2),rep("labeled",n3)))
set.seed(2017)
folds = c(rep(1:10,n/10+1))
folds = folds[-1]
folds = folds[sample(1:length(folds))]
tb = 0
for(kk in 1:10){
  ind_test = folds == kk
  train = dat2[-ind_test,]
  test = dat2[ind_test,]
  rf.fit = randomForest(cl~.,data = train)
  pred = predict(rf.fit, newdata = test, type = "class")
  tb = tb + table(pred,test$cl)
}
acc0 = sum(diag(tb))/sum(tb)
f1score0 = 2*tb[1,1]/(tb[1,2]+tb[2,1]+2*tb[1,1])

alpha_est_tcdb = matrix(0, nrow=3,ncol=5)
alpha_est_tcdb[1,] = c(round(alpha,3), round(rst3$alpha,3), round(rst2$alpha,3), round(rst4$alpha,3),round(rst5$alpha,3))
alpha_est_tcdb[2,] = c(round(acc0, 3),round(acc3, 3),round(acc2, 3),round(acc4, 3),round(acc5, 3)) 
alpha_est_tcdb[3,] = c(round(f1score0,3),round(f1score3,3),round(f1score2,3),round(f1score4,3),round(f1score5,3))
colnames(alpha_est_tcdb) <- c("ideal", "C-patra/sen","C-roc","ROC","SPY")
rownames(alpha_est_tcdb) <- c('alpha','accuracy','F1 score')

## print result to a tex file
print(xtable(alpha_est_tcdb, caption = "Comparison of methods for protein signaling data.", 
             label = "tab:tcdb_swissprot"),type="latex",
      file="figs/comparison_varying_data.tex", include.rownames = T)

#file.copy(from=list.files('figs',full.names=TRUE),to='../ms/figs/',overwrite=TRUE)




