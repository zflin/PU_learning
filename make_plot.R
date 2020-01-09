if(!dir.exists("figs")){
  dir.create("figs")
}

########################################################################################
#### varying feature 
########################################################################################
rm(list=ls())

load('results/waveform_vary_feature.RData')

pdf("figs/comparison_varying_feature.pdf", width = 12, height = 6)

layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE), 
       widths=c(2,1), heights=c(2,2))
par(mar=c(5,6,1,5))
plot(1:nfeat,alphas_all,
     ylim=c(0,1),pch=4,col=2,
     xlab='Feature',ylab="",cex.lab=2,cex.axis=1.9)
abline(h=alpha, col=4, lty=3)
par(new = T)
plot(1:nfeat, imp, type = 'h', lty=2, axes=F, xlab=NA, ylab=NA)
axis(side=4,cex.axis=1.9)
mtext(side=4,line=3,'Importance',cex=2)
mtext(side=2,line=3,expression(widehat(alpha)[0]),cex=2)

legend("topright", legend = c('True Alpha','Importance','PS'),
       pch = c(NA, NA, 4), lty = c(3, 2, NA),
       col=c(4,1,2),cex=2)
cex_meth <- 1.8
st <- -2.2
sp <- 1.9
mtext(substitute(paste("True alpha = ", m), list(m=round(alpha,3))), 3, line = st, adj=0.01, cex=cex_meth)
mtext(substitute(paste("C-PS = ", m), list(m=round(alpha3,3))), 3, line = st-sp, adj=0.01, cex=cex_meth)
mtext(substitute(paste("C-ROC = ", m), list(m=round(alpha2,3))), 3, line = st-2*sp, adj=0.01, cex=cex_meth)
mtext(substitute(paste("ROC = ", m), list(m=round(alpha4,3))), 3, line = st-3*sp, adj=0.01, cex=cex_meth)
mtext(substitute(paste("SPY = ", m), list(m=round(alpha5,3))), 3, line = st-4*sp, adj=0.01, cex=cex_meth)
mtext(substitute(paste("max(PS) = ", m), list(m=round(max(alphas_all),3))), 3, line = st-5*sp, adj=0.01, cex=cex_meth)

par(mar=c(2.5,1,4,1))

lwd <- 4

kde_cex <- 1.9
plot(density(tmp0_5), main = "KDEs of 5th feature", axes=F, xlab=NA, ylab=NA, ylim=c(0,0.3),cex.main=kde_cex,lwd=lwd)
lines(density(tmp1_5),col=2,lwd=lwd,lty=2)
# axis(side=4)
# mtext(side = 4, line = 3, 'Density')
legend("topleft",c("U","L"),lty=c(1,2), col=1:2, cex=2, bty = "n",lwd=lwd)

plot(density(tmp0_8), ylim=c(0,0.35), main = "KDEs of 8th feature", axes=F, xlab=NA, ylab=NA,cex.main=kde_cex,lwd=lwd)
lines(density(tmp1_8),col=2,lwd=lwd,lty=2)
# axis(side=4)
# mtext(side = 4, line = 3, 'Density')
legend("topright",c("U","L"),lty=c(1,2), col=1:2, cex=2,bty = "n",lwd=lwd)
dev.off()

########################################################################################
#### varying alpha 
########################################################################################
rm(list=ls())

load('results/waveform_vary_alpha_croc_patrasen.RData')
load('results/waveform_vary_alpha_spy.RData')
load('results/waveform_vary_alpha_roc.RData')

alphas = alpha_est_roc$alpha

pdf("figs/comparison_varying_alpha.pdf",width = 10, height = 3.5)

cols <- c("#000000", "#E69F00", "#56B4E9", "#009E73")

par(mfrow=c(1,3))
plot(alphas, alpha_est_croc_patrasen$alpha2, ylim = c(0,1), pch=19, cex=0.85,
     ylab = expression(widehat(alpha)[0]), xlab = expression(alpha),col=cols[1])
points(alphas, alpha_est_croc_patrasen$alpha3, col=cols[2], pch=2)
points(alphas, alpha_est_roc$alpha4, col=cols[3], pch=3)
points(alphas, alpha_est_spy$alpha5, col=cols[4], pch=4)
abline(0,1, lty=2)
legend('bottomright', legend = c('C-PS', 'C-ROC', 'ROC', 'SPY'), pch=c(2,19,3,4), col=c(cols[2],cols[1],cols[3],cols[4]))

plot(alphas, alpha_est_croc_patrasen$acc2, ylim = c(0,1), pch=19, cex=0.85,
     ylab = 'accuracy', xlab = expression(alpha),col=cols[1])
points(alphas, alpha_est_croc_patrasen$acc3, col=cols[2], pch=2)
points(alphas, alpha_est_roc$acc4, col=cols[3], pch=3)
points(alphas, alpha_est_spy$acc5, col=cols[4], pch=4)
legend('bottomright', legend = c('C-PS', 'C-ROC', 'ROC', 'SPY'), pch=c(2,19,3,4), col=c(cols[2],cols[1],cols[3],cols[4]))

plot(alphas, alpha_est_croc_patrasen$f1score2, ylim = c(0,1), pch=19, cex=0.85,
     ylab = 'F1 score', xlab = expression(alpha),col=cols[1])
points(alphas, alpha_est_croc_patrasen$f1score3, col=cols[2], pch=2)
points(alphas, alpha_est_roc$f1score4, col=cols[3], pch=3)
points(alphas, alpha_est_spy$f1score5, col=cols[4], pch=4)
legend('bottomright', legend = c('C-PS', 'C-ROC', 'ROC', 'SPY'), pch=c(2,19,3,4), col=c(cols[2],cols[1],cols[3],cols[4]))
par(mfrow=c(1,1))

dev.off()


########################################################################################
#### varying sample size 
########################################################################################
rm(list=ls())

load('results/waveform_vary_size_croc_patrasen.RData')
load('results/waveform_vary_size_roc.RData')
load('results/waveform_vary_size_spy.RData')

pdf("figs/comparison_varying_samplesize.pdf", height = 6, width = 8)

alphas = c(.1, .5, .9)
par(mfrow=c(3,4))
opar = par( oma=c(2,3.7,2,0), mar=c(1.5,1,1,1.5)+0.1 )
for(ii in 1:length(alphas)){
  tmp = alpha_est_croc_patrasen[alpha_est_croc_patrasen$alpha==alphas[ii],]
  boxplot(tmp$alpha3~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F, yaxt="n")
  axis(side=2, at=c(alphas), labels=c(alphas))
  abline(h=alphas[ii], col=2)
  if(ii==1){mtext( 'C-PS', side=3, line=0.25, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=1, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
  
  tmp = alpha_est_croc_patrasen[alpha_est_croc_patrasen$alpha==alphas[ii],]
  boxplot(tmp$alpha2~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F, yaxt="n")
  axis(side=2, at=c(alphas), labels=c(alphas))
  abline(h=alphas[ii], col=2)
  if(ii==2){mtext( expression(widehat(alpha)[0]), side=2, line=2, at=grconvertY(0.5,'npc','nic'), outer=TRUE )}
  if(ii==1){mtext( 'C-ROC', side=3, line=0.25, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=1, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
  
  tmp = alpha_est_roc[alpha_est_roc$alpha==alphas[ii],]
  boxplot(tmp$alpha4~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F, yaxt="n")
  axis(side=2, at=c(alphas), labels=c(alphas))
  abline(h=alphas[ii], col=2)
  if(ii==1){mtext( 'ROC', side=3, line=0.25, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=1, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
  
  tmp = alpha_est_spy[alpha_est_spy$alpha==alphas[ii],]
  boxplot(tmp$alpha5~tmp$grp, ylim=c(0,1), names=names(table(tmp$n0)), ann=F, yaxt="n")
  axis(side=2, at=c(alphas), labels=c(alphas))
  abline(h=alphas[ii], col=2)
  if(ii==1){mtext( 'SPY', side=3, line=0.25, at=grconvertX(0.5,'npc','nic'), outer=TRUE, font=2 )}
  if(ii==3){mtext( 'n', side=1, line=1, at=grconvertX(0.5,'npc','nic'), outer=TRUE )}
}
par(opar)

dev.off()

#####

