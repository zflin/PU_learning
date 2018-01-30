rm(list=ls())
set.seed(1234)
library(randomForest)
library(sm) ## kde 
source("funcs.R")
source("funcs_patrasen.R")
source("funcs_scott.R")
library(xtable)

#### varying feature ######

# add_legend <- function(...) {
#   opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
#               mar=c(0, 0, 0, 0), new=TRUE)
#   on.exit(par(opar))
#   plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
#   legend(...)
# }

load('data/waveform_vary_feature.RData')

pdf("figs/comparison_varying_feature.pdf", width = 10, height = 6)
layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE), 
       widths=c(2,1), heights=c(2,2))
plot(1:nfeat, alphas_all, ylim = c(0,1), pch=4, col=2, xlab = 'feature', ylab = expression(hat(alpha)[0]))
points(1:nfeat, imp/max(imp), type = 'h', lty=2)
# points(1:nfeat, alphas_all[,2],col=2,pch=2)
# points(1:nfeat, alphas_all[,3],col=4,pch=4)
abline(h=alpha, col=4, lty=3)
# abline(h=alpha2, col=2, lty=2)
# abline(h=alpha3, col=3, lty=3)
# abline(h=alpha4, col=4, lty=4)
# legend('topleft', legend = c('C-storey','C-patra/sen','patra/sen'),
#        col=c(1,2,4),pch=c(1,2,4))
legend("topright", legend = c('true alpha','importance','patra/sen'),
       pch = c(NA, NA, 4), lty = c(3, 2, NA),
       col=c(4,1,2))
mtext(substitute(paste("true alpha = ", m), list(m=round(alpha,3))), 3, line = -2, adj=0.01, cex=0.9)
mtext(substitute(paste("C-storey = ", m), list(m=round(alpha2,3))), 3, line = -3, adj=0.01, cex=0.9)
mtext(substitute(paste("C-patra/sen = ", m), list(m=round(alpha3,3))), 3, line = -4, adj=0.01, cex=0.9)
mtext(substitute(paste("scott = ", m), list(m=round(alpha4,3))), 3, line = -5, adj=0.01, cex=0.9)
mtext(substitute(paste("spy = ", m), list(m=round(alpha5,3))), 3, line = -6, adj=0.01, cex=0.9)
mtext(substitute(paste("max(patra/sen) = ", m), list(m=round(max(alphas_all),3))), 3, line = -7, adj=0.01, cex=0.9)

plot(density(tmp0_5), main = "KDEs of 5th feature", xlab = "", ylim=c(0,0.3))
lines(density(tmp1_5),col=2)
legend("topright",c("unlabeled","labeled"),lty=c(1,1), col=1:2, cex=0.7, bty = "n")
# add_legend("top", legend=c("unlabeled","labeled"),lty=c(1,1), col=1:2,
#            horiz=TRUE, bty='n', cex=1)

plot(density(tmp0_8), ylim=c(0,0.35), main = "KDEs of 8th feature", xlab = "")
lines(density(tmp1_8),col=2)
legend("topright",c("unlabeled","labeled"),lty=c(1,1), col=1:2, cex=0.7,bty = "n")

# patra16 = EstimateParams(tmp0, ecdf(tmp1),  alpha_grid = (1:1000)/1000)
# dder = ComputeSecondDeriv(patra16$out$distance)
# plot(patra16$out$alpha_grid,patra16$out$distance, xlab = 'alpha0', ylab = 'distance',
#      main = "patra/sen diagnostic plot of 16th feature")
# lines(patra16$out$alpha_grid ,dder*(max(patra16$out$distance)/max(dder)),col='red')
# legend("top",c("Distance","Scaled 2nd derivative"),
#        lty=c(1,1), col=c("blue","red"), cex=1)

dev.off()

#### varying alpha #########

load('data/waveform_vary_alpha_storey_patrasen_spy.RData')
load('data/waveform_vary_alpha_scott.RData')

# alpha_est$acc2[alpha_est$acc2<0.5] <- 1-alpha_est$acc2[alpha_est$acc2<0.5]
# alpha_est$acc3[alpha_est$acc3<0.5] <- 1-alpha_est$acc3[alpha_est$acc3<0.5]
# alpha_est$acc4[alpha_est$acc4<0.5] <- 1-alpha_est$acc4[alpha_est$acc4<0.5]

alphas = alpha_est_scott$alpha

pdf("figs/comparison_varying_alpha.pdf",width = 10, height = 3.5)
par(mfrow=c(1,3))
plot(alphas, alpha_est_storey_patrasen_spy$alpha2, ylim = c(0,1), 
     ylab = expression(hat(alpha)[0]), xlab = expression(alpha))
points(alphas, alpha_est_storey_patrasen_spy$alpha3, col=2, pch=2)
points(alphas, alpha_est_scott$alpha4, col=4, pch=3)
points(alphas, alpha_est_storey_patrasen_spy$alpha5, col=6, pch=4)
abline(0,1)
legend('bottomright', legend = c('C-storey', 'C-patra/sen', 'scott', 'spy'), pch=1:4, col=c(1,2,4,6))

plot(alphas, alpha_est_storey_patrasen_spy$acc2, ylim = c(0,1), ylab = 'accuracy', xlab = expression(alpha))
points(alphas, alpha_est_storey_patrasen_spy$acc3, col=2, pch=2)
points(alphas, alpha_est_scott$acc4, col=4, pch=3)
points(alphas, alpha_est_storey_patrasen_spy$acc5, col=6, pch=4)
legend('bottomright', legend = c('C-storey', 'C-patra/sen', 'scott', 'spy'), pch=1:4, col=c(1,2,4,6))

plot(alphas, alpha_est_storey_patrasen_spy$f1score2, ylim = c(0,1), ylab = 'F1 score', xlab = expression(alpha))
points(alphas, alpha_est_storey_patrasen_spy$f1score3, col=2, pch=2)
points(alphas, alpha_est_scott$f1score4, col=4, pch=3)
points(alphas, alpha_est_storey_patrasen_spy$f1score5, col=6, pch=4)
legend('bottomright', legend = c('C-storey', 'C-patra/sen', 'scott', 'spy'), pch=1:4, col=c(1,2,4,6))
par(mfrow=c(1,1))
dev.off()

### varying sample size #####

load('data/waveform_vary_size_patra_storey.RData')
load('data/waveform_vary_size_scott.RData')
load('data/waveform_vary_size_spy.RData')

pdf("figs/comparison_varying_samplesize.pdf")
alphas = c(.1, .5, .9)
par(mfrow=c(3,4))
par( oma=c(6,6,6,0), mar=c(2,2,1,1)+0.1 )
for(ii in 1:length(alphas)){
  tmp = alpha_est_patra_storey[alpha_est_patra_storey$alpha==alphas[ii],]
  boxplot(tmp$alpha2~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F)
  abline(h=alphas[ii], col=2)
  if(ii==2){mtext( 'alpha0 estimation', side=2, line=2, at=grconvertY(0.5,'npc','nic'), outer=TRUE )}
  if(ii==1){mtext( 'C-storey', side=3, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
  
  tmp = alpha_est_patra_storey[alpha_est_patra_storey$alpha==alphas[ii],]
  boxplot(tmp$alpha3~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F)
  abline(h=alphas[ii], col=2)
  if(ii==1){mtext( 'C-patra/sen', side=3, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
  
  tmp = alpha_est_scott[alpha_est_scott$alpha==alphas[ii],]
  boxplot(tmp$alpha4~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F)
  abline(h=alphas[ii], col=2)
  if(ii==1){mtext( 'scott', side=3, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
  
  tmp = alpha_est_spy[alpha_est_spy$alpha==alphas[ii],]
  boxplot(tmp$alpha5~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F)
  abline(h=alphas[ii], col=2)
  if(ii==1){mtext( 'spy', side=3, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
}
dev.off()

#####

file.copy(from=list.files('figs',full.names=TRUE),to='../ms/figs/',overwrite=TRUE)


