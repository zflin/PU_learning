if(!dir.exists("figs")){
  dir.create("figs")
}

########################################################################################
#### varying feature 
########################################################################################
rm(list=ls())

load('data/waveform_vary_feature.RData')

pdf("figs/comparison_varying_feature.pdf", width = 10, height = 6)

layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE), 
       widths=c(2,1), heights=c(2,2))
plot(1:nfeat, alphas_all, ylim = c(0,1), pch=4, col=2, xlab = 'feature', ylab = expression(widehat(alpha)[0]))
abline(h=alpha, col=4, lty=3)
par(new = T)
plot(1:nfeat, imp, type = 'h', lty=2, axes=F, xlab=NA, ylab=NA)
axis(side=4)
mtext(side = 4, line = 3, 'importance')

legend("topright", legend = c('true alpha','importance','patra/sen'),
       pch = c(NA, NA, 4), lty = c(3, 2, NA),
       col=c(4,1,2))
mtext(substitute(paste("true alpha = ", m), list(m=round(alpha,3))), 3, line = -2, adj=0.01, cex=0.9)
mtext(substitute(paste("C-patra/sen = ", m), list(m=round(alpha3,3))), 3, line = -3, adj=0.01, cex=0.9)
mtext(substitute(paste("C-roc = ", m), list(m=round(alpha2,3))), 3, line = -4, adj=0.01, cex=0.9)
mtext(substitute(paste("ROC = ", m), list(m=round(alpha4,3))), 3, line = -5, adj=0.01, cex=0.9)
mtext(substitute(paste("SPY = ", m), list(m=round(alpha5,3))), 3, line = -6, adj=0.01, cex=0.9)
mtext(substitute(paste("max(patra/sen) = ", m), list(m=round(max(alphas_all),3))), 3, line = -7, adj=0.01, cex=0.9)

plot(density(tmp0_5), main = "KDEs of 5th feature", axes=F, xlab=NA, ylab=NA, ylim=c(0,0.3))
lines(density(tmp1_5),col=2)
# axis(side=4)
# mtext(side = 4, line = 3, 'Density')
legend("topright",c("unlabeled","labeled"),lty=c(1,1), col=1:2, cex=0.7, bty = "n")

plot(density(tmp0_8), ylim=c(0,0.35), main = "KDEs of 8th feature", axes=F, xlab=NA, ylab=NA)
lines(density(tmp1_8),col=2)
# axis(side=4)
# mtext(side = 4, line = 3, 'Density')
legend("topright",c("unlabeled","labeled"),lty=c(1,1), col=1:2, cex=0.7,bty = "n")

dev.off()

########################################################################################
#### varying alpha 
########################################################################################
rm(list=ls())

load('data/waveform_vary_alpha_croc_patrasen.RData')
load('data/waveform_vary_alpha_spy.RData')
load('data/waveform_vary_alpha_roc.RData')

alphas = alpha_est_roc$alpha

pdf("figs/comparison_varying_alpha.pdf",width = 10, height = 3.5)

par(mfrow=c(1,3))
plot(alphas, alpha_est_croc_patrasen$alpha2, ylim = c(0,1), pch=19, cex=0.7,
     ylab = expression(widehat(alpha)[0]), xlab = expression(alpha))
points(alphas, alpha_est_croc_patrasen$alpha3, col=2, pch=2)
points(alphas, alpha_est_roc$alpha4, col=4, pch=3)
points(alphas, alpha_est_spy$alpha5, col=6, pch=4)
abline(0,1, lty=2)
legend('bottomright', legend = c('C-patra/sen', 'C-roc', 'ROC', 'SPY'), pch=c(2,19,3,4), col=c(2,1,4,6))

plot(alphas, alpha_est_croc_patrasen$acc2, ylim = c(0,1), pch=19, cex=1,
     ylab = 'accuracy', xlab = expression(alpha))
points(alphas, alpha_est_croc_patrasen$acc3, col="#FF0000C0", pch=2)
points(alphas, alpha_est_roc$acc4, col=4, pch=3)
points(alphas, alpha_est_spy$acc5, col=6, pch=4)
legend('bottomright', legend = c('C-patra/sen', 'C-roc', 'ROC', 'SPY'), pch=c(2,19,3,4), col=c(2,1,4,6))

plot(alphas, alpha_est_croc_patrasen$f1score2, ylim = c(0,1), pch=19, cex=1,
     ylab = 'F1 score', xlab = expression(alpha))
points(alphas, alpha_est_croc_patrasen$f1score3, col="#FF0000C0", pch=2)
points(alphas, alpha_est_roc$f1score4, col=4, pch=3)
points(alphas, alpha_est_spy$f1score5, col=6, pch=4)
legend('bottomright', legend = c('C-patra/sen', 'C-roc', 'ROC', 'SPY'), pch=c(2,19,3,4), col=c(2,1,4,6))
par(mfrow=c(1,1))

dev.off()


########################################################################################
#### varying sample size 
########################################################################################
rm(list=ls())

load('data/waveform_vary_size_croc_patrasen.RData')
load('data/waveform_vary_size_roc.RData')
load('data/waveform_vary_size_spy.RData')

pdf("figs/comparison_varying_samplesize.pdf", height = 6, width = 8)

alphas = c(.1, .5, .9)
par(mfrow=c(3,4))
opar = par( oma=c(3,3.7,5,0), mar=c(1.5,1,1,1.5)+0.1 )
for(ii in 1:length(alphas)){
  tmp = alpha_est_croc_patrasen[alpha_est_croc_patrasen$alpha==alphas[ii],]
  boxplot(tmp$alpha3~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F, yaxt="n")
  axis(side=2, at=c(alphas), labels=c(alphas))
  abline(h=alphas[ii], col=2)
  if(ii==1){mtext( 'C-patra/sen', side=3, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
  
  tmp = alpha_est_croc_patrasen[alpha_est_croc_patrasen$alpha==alphas[ii],]
  boxplot(tmp$alpha2~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F, yaxt="n")
  axis(side=2, at=c(alphas), labels=c(alphas))
  abline(h=alphas[ii], col=2)
  if(ii==2){mtext( expression(widehat(alpha)[0]), side=2, line=2, at=grconvertY(0.5,'npc','nic'), outer=TRUE )}
  if(ii==1){mtext( 'C-roc', side=3, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
  
  tmp = alpha_est_roc[alpha_est_roc$alpha==alphas[ii],]
  boxplot(tmp$alpha4~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F, yaxt="n")
  axis(side=2, at=c(alphas), labels=c(alphas))
  abline(h=alphas[ii], col=2)
  if(ii==1){mtext( 'ROC', side=3, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
  
  tmp = alpha_est_spy[alpha_est_spy$alpha==alphas[ii],]
  boxplot(tmp$alpha5~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F, yaxt="n")
  axis(side=2, at=c(alphas), labels=c(alphas))
  abline(h=alphas[ii], col=2)
  if(ii==1){mtext( 'SPY', side=3, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=2, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
}
par(opar)

dev.off()

#####

# file.copy(from=list.files('figs',full.names=TRUE),to='../ms/figs/',overwrite=TRUE)


