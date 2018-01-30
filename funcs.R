
## storey's estimation
##' p0: probability prediction of unlabeled samples
##' p1: probability prediction of positive samples
##' cl: true class label of unlabeled samples, default is NULL, only used to get accuracy
##' alpha: true alpha, default is 0.5, only used to plot
##' plotit: to plot or not
# method_c_storey <- function(p0, p1, cl=NULL, alpha=0.5, plotit=F){
#   pi_hat = length(p1)/(length(p0)+length(p1)) # proportion of positive samples
#   pi_prime = (1-pi_hat)/pi_hat
#   tstar = seq(min(p0),max(p0), length.out = length(unique(p0))) # resample probability predictions
#   #tstar = sort(p0)
#   G0 = ecdf(p0) # cdf of unlabeled prediction
#   G1 = ecdf(p1) # cdf of labeled prediction
#   G0d = G0(tstar) # realization of ecdf
#   G1d = G1(tstar) #
#   
#   kt = (G0d-G1d)/(1-G1d)
#   V = (G0d - (1+pi_prime)*G0d*G1d + pi_prime*G1d)*(1-G0d)/(1-G1d)^3/length(p0) # variance
#   
#   ## to remove NA in kt (and the related tstar and V)
#   ind = is.na(kt)
#   kt = kt[!ind]
#   tstar = tstar[!ind]
#   V = V[!ind]
#   ## to remove the first negative value and the values after in kt (and the related tstar and V)
#   ind = cumsum(kt < 0)
#   if(min(ind)==1){
#     ind1 = length(ind)
#   }else{
#     ind1 = max(which(ind==0))
#   }
#   kt = kt[1:(ind1-1)]
#   tstar = tstar[1:(ind1-1)]
#   V = V[1:(ind1-1)]
#   G0d = G0(tstar)
#   G1d = G1(tstar)
#   
#   ## 95% lower fonfidence band
#   lt = kt - qnorm(0.95)*sqrt(V)
#   
#   ## kde of p1 and p0
#   g1 = sm.density(p1, eval.points=tstar, display = "none")
#   g0 = sm.density(p0, eval.points=tstar, display = "none")
#   
#   ## 
#   dt = 1-g0$estimate/g1$estimate
#   
#   dt1 = (G0d-G1d)*(g1$estimate)
#   dt2 = (g1$estimate-g0$estimate)*(1-G1d)
#   Dt =  (dt1 + dt2)/2
#   #Dt = (Dt-min(Dt)+1)/(max(Dt)-min(Dt)+1) ## scale for good-looking plot
#   
#   alpha2 = kt[which.max(Dt)] ## alpha0 estimate
#   
#   if(plotit){
#     ## empirical cdf of p1 and p0
#     plot(G0, main="ecdf")
#     lines(G1,col=2)
#     legend('bottomright', legend = c('unlabeled', 'labeled'), col=1:2, lty = 1)
#     
#     ## numerator and denominator of kt
#     plot(tstar,G0d-G1d, main="ecdf", ylim = c(-0.5,1))
#     lines(tstar,1-G1d,col=2)
#     abline(h=0,lty=2,col=4)
#     legend('topright', legend = c('unlabeled', 'labeled'), col=1:2, lty = 1)
#     
#     ## 
#     opar = par(mar=c(5.1, 4.1, 4.1, 8.1))
#     plot(tstar, dt1, ylim = c(0,1), ylab = "", xlab = "t", main = "c-storey")
#     points(tstar, dt2, col=4)
#     points(tstar, Dt, col=3)
#     points(tstar,kt,col=2)
#     abline(h=alpha)
#     abline(h=0)
#     abline(v = tstar[which.max((kt+dt)*g1$estimate*(1-G1d)/2)], lty=2)
#     legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, inset=c(-0.25,0), legend=c("k(t)","d_1(t)","d_2(t)", "D(t)"), 
#            pch=c(19,19,19,19), col=c(2,1,4,3))
#     par(opar)
#     
#     # colvec = c("#00000060","#99000060",8,4,6)
#     # opar = par(mar=c(5.1, 4.1, 4.1, 8.1))
#     # plot(tstar,kt, ylim = c(0,1), pch=19, cex=0.5, col=colvec[1],
#     #      ylab = "k(t)", xlab = "t")
#     # points(tstar,lt, col=colvec[2], pch=19, cex=0.3)
#     # points(tstar,dt, col=colvec[3], pch=19, cex=0.3)
#     # lines(tstar,dg, col=colvec[4])
#     # abline(h=alpha, col=colvec[5], xlim=c(0,max(tstar)))
#     # abline(v=tstar[which.min(dg)], lty=2)
#     # legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, inset=c(-0.25,0), legend=c("k(t)","l(t)","d(t)","D(t)","true alpha"), 
#     #        pch=c(19,19,19,NA,NA,NA), lty=c(NA,NA,NA,1,1), cex=0.7,
#     #        col=colvec)
#     # par(opar)
#   }
#   
#   threshold = pi_hat/(pi_hat + 2*(1-pi_hat)*(1-alpha2))
#   tb2 = table(p0 > threshold,cl) # confusion matrix
#   acc2 = sum(diag(tb2))/sum(tb2)
#   return(list(alpha=alpha2,acc=acc2))
# }

## storey's estimation
##' p0: probability prediction of unlabeled samples
##' p1: probability prediction of positive samples
##' cl: true class label of unlabeled samples, default is NULL, only used to get accuracy
##' alpha: true alpha, default is 0.5, only used to plot
##' plotit: to plot or not
method_c_storey <- function(p0, p1, cl=NULL, alpha=0.5, plotit=F){
  pi_hat = length(p1)/(length(p0)+length(p1)) # proportion of positive samples
  pi_prime = (1-pi_hat)/pi_hat
  tstar = seq(min(p0),max(p0), length.out = 500) # resample probability predictions
  #tstar = sort(p0)
  G0 = ecdf(p0) # cdf of unlabeled prediction
  G1 = ecdf(p1) # cdf of labeled prediction
  G0d = G0(tstar) # realization of ecdf
  G1d = G1(tstar) #
  
  kt = (G0d-G1d)/(1-G1d)
  V = (G0d - (1+pi_prime)*G0d*G1d + pi_prime*G1d)*(1-G0d)/(1-G1d)^3/length(p0) # variance
  
  ## to remove NA in kt (and the related tstar and V)
  ind = is.na(kt)
  kt = kt[!ind]
  tstar = tstar[!ind]
  V = V[!ind]
  ## to remove the first negative value and the values after in kt (and the related tstar and V)
  ind = cumsum(kt < 0)
  if(min(ind)==1){
    ind1 = length(ind)
  }else{
    ind1 = max(which(ind==0))
  }
  kt = kt[1:(ind1-1)]
  tstar = tstar[1:(ind1-1)]
  V = V[1:(ind1-1)]
  G0d = G0(tstar)
  G1d = G1(tstar)
  
  # ## 95% lower fonfidence band
  # lt = kt - qnorm(0.95)*sqrt(V)
  # 
  # ## kde of p1 and p0
  # g1 = sm.density(p1, eval.points=tstar, display = "none")
  # g0 = sm.density(p0, eval.points=tstar, display = "none")
  # 
  # ## 
  # dt = 1-g0$estimate/g1$estimate
  # 
  # ## Dt idea
  # dt1 = (G0d-G1d)*(g1$estimate)
  # dt2 = (g1$estimate-g0$estimate)*(1-G1d)
  # Dt =  (dt1 + dt2)/2
  # alpha.Dt = kt[which.max(Dt)] ## alpha0 estimate
  
  ## second derivative idea
  df0 = 2
  y =  diff(G0d-G1d, differences = df0)/diff(1-G1d, differences = df0)
  y[!is.finite(y)] = 0
  y[y >= 1-0.001] = 0
  alpha0 = kt[which.max(y)]
  
  if(plotit){
    plot(tstar,kt, ylim=c(0,1), cex=0.1)
    abline(h=alpha)
    #abline(v = tstar[which.max(Dt)], lty=2)
    #lines(tstar,dt)
    points(tstar[(df0+1):length(tstar)], y,col=2)
    abline(v=tstar[which.max(y)])
  }
  
  threshold = pi_hat/(pi_hat + 2*(1-pi_hat)*(1-alpha0))
  decision = p0 > threshold
  tb = table(decision,cl) # confusion matrix
  tb = as.vector(tb)
  if(length(unique(decision))==1){
    if(unique(decision)){
      tb = c(0,0,tb)
    }else{
      tb = c(tb,0,0)
    }
  }
  acc = sum(tb[c(1,4)])/sum(tb)
  #f1score = 2*tb[2,2]/(tb[1,2]+tb[2,1]+2*tb[2,2])
  return(list(alpha=alpha0,acc=acc, tb=tb))
}

## patra/sen's estimation
##' p0: probability prediction of unlabeled samples
##' p1: probability prediction of positive samples
##' cl: true class label of unlabeled samples, default is NULL, only used to get accuracy
##' alpha: true alpha, default is 0.5, only used to plot
##' plotit: to plot or not
method_c_patrasen <- function(p0, p1, cl=NULL, alpha=0.5, plotit=TRUE){
  pi_hat = length(p1)/(length(p0)+length(p1))
  out <- EstMixMdl(p0, ecdf(p1),alpha_grid = (1:1000)/1000)
 
  dder = ComputeSecondDeriv(out$distance)
  ix <- which.max(dder)
  alpha0 <- out$alpha_grid[ix]
  
  # tmp = supsmu(out$alpha_grid,out$distance)
  # dder = ComputeSecondDeriv(tmp$y)
  # ix <- which.max(dder)
  # alpha3 <- tmp$x[ix]
  
  if(plotit){
    # plot(tmp, xlab = 'alpha0', ylab = 'distance')
    # lines(tmp$x ,dder*(max(tmp$y)/max(dder)),col='red')
    # legend("topright",c("Distance","Scaled 2nd derivative"),
    #        lty=c(1,1), col=c("blue","red"), cex=1)
    plot(out$alpha_grid,out$distance, xlab = 'alpha0', ylab = 'distance')
    lines(out$alpha_grid ,dder*(max(out$distance)/max(dder)),col='red')
    legend("topright",c("Distance","Scaled 2nd derivative"),
           lty=c(1,1), col=c("blue","red"), cex=1)
  }
  
  threshold = pi_hat/(pi_hat + 2*(1-pi_hat)*(1-alpha0))
  decision = p0 > threshold
  tb = table(decision,cl) # confusion matrix
  tb = as.vector(tb)
  if(length(unique(decision))==1){
    if(unique(decision)){
      tb = c(0,0,tb)
    }else{
      tb = c(tb,0,0)
    }
  }
  acc = sum(tb[c(1,4)])/sum(tb)
  #f1score = 2*tb[2,2]/(tb[1,2]+tb[2,1]+2*tb[2,2])
  return(list(alpha=alpha0,acc=acc, tb=tb))
}

## scott's estimation
##' nrep = 5: number of repeatitions to get robust result
method_scott <- function(x0, x1, cl=NULL, alpha=0.5, nrep = 5){
  pi_hat = nrow(x1)/(nrow(x0)+nrow(x1))
  pi1t_roc = rep(0,nrep)
  pi0t_roc = rep(0,nrep)
  yhat = rep(0,nrow(x0))
  for(ii in 1:nrep){
    cpe.fit = cpe(x1, x0)
    yhat = yhat + cpe.fit$Y_fit$probabilities[1:nrow(x0),1]
    eta = cpe.fit$eta
    y = cpe.fit$labels
    nvec = cpe.fit$nvec
    pi1t_roc[ii] = kappahat(eta, y, nvec[4], nvec[3])
    pi0t_roc[ii] = kappahat(1-eta, -y, nvec[3], nvec[4])
  }
  yhat = yhat/nrep
  
  mpi1t_roc = median(pi1t_roc)
  mpi0t_roc = median(pi0t_roc)
  pi1_roc = 1 - mpi1t_roc*(1-mpi0t_roc)/(1-mpi0t_roc*mpi1t_roc)
  alpha0 = 1 - mpi0t_roc*(1-mpi1t_roc)/(1-mpi0t_roc*mpi1t_roc)
  
  decision = yhat > pi_hat/(pi_hat + 2*(1-pi_hat)*(1-alpha0))
  tb = table(decision,cl) # confusion matrix
  tb = as.vector(tb)
  if(length(unique(decision))==1){
    if(unique(decision)){
      tb = c(0,0,tb)
    }else{
      tb = c(tb,0,0)
    }
  }
  acc = sum(tb[c(1,4)])/sum(tb)
  #f1score = 2*tb[2,2]/(tb[1,2]+tb[2,1]+2*tb[2,2])
  return(list(alpha=alpha0,acc=acc, tb=tb))
}

## spy tech
##' s = proportion of positive as spy set
##' l = noise level to determine threshold for reliable negative set
method_spy <- function(x0,x1, cl=NULL, s=0.15,l=0.15){
  indspy <- ceiling(nrow(x1)*(1-s)):nrow(x1)
  datspy = rbind(x1,x0)
  datspy$cl <- as.factor(c(rep("labeled",nrow(x1)), rep("unlabeled",nrow(x0))))
  datspy$cl[indspy] = "unlabeled"
  
  set.seed(17)
  rf.fit1 <- randomForest(cl~.,data=datspy)
  rf.preds1 <- predict(rf.fit1,type='prob')
  preds1 <- rf.preds1[,1]
  lu1 <- rep(0,length(preds1))
  lu1[datspy$cl=="labeled"] <- 1
  p01 = preds1[lu1==0]
  p11 = preds1[lu1==1]
  
  pspy = preds1[indspy]
  
  thres = min(quantile(pspy, l))
  xrn = x0[preds1[(nrow(x0)+1):length(preds1)]<thres,] # reliable negative
  dat_PRN = rbind(x1, xrn) # P + RN
  dat_PRN$cl <- as.factor(c(rep("labeled",nrow(x1)), rep("unlabeled",nrow(xrn))))
  
  #set.seed(17)
  rf.fit2 <- randomForest(cl~.,data=dat_PRN)
  rf.preds2 <- predict(rf.fit2, x0)
  tb = table(rf.preds2,1-cl)
  
  alpha0 = mean(rf.preds2=="unlabeled")
  acc = sum(diag(tb))/sum(tb)
  #f1score = 2*tb[2,2]/(tb[1,2]+tb[2,1]+2*tb[2,2])
  return(list(alpha=alpha0,acc=acc, tb=as.vector(tb)))
}




