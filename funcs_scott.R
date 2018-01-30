library(MASS)
library(pdist)
library(LiblineaR)

transformFeatures <- function(X, Omega, beta){
  ##' X: n x p 
  ##' Omega: p x D 
  ##' beta: 1 x D
  ##' 
  ##' Z: n x D
  D = ncol(Omega)
  Z = cos(sweep(X %*% Omega,2,beta, '+'))*sqrt(2/D)
  return(Z)
}

auc <- function(eta, Y, res=20){
  ##' eta: vector of probabilities of a positive
  ##' Y: vector of true labels
  ##' res: number of threshold points to test, defaults to 20
  ##' 
  ##' area: AUC
  ##' alpha: vector of probabilities of a false positive
  ##' beta: vector of probabilities of a true positive
  n_pos = sum(Y==1)
  n_neg = sum(Y==-1)
  alpha = rep(0,res)
  beta = rep(0,res)
  thresh = seq(1,0,length.out = res)
  for(ii in 1:res){
    t = thresh[ii]
    alpha[ii] = sum((eta>t)*(Y==-1))/n_neg
    beta[ii] = sum((eta>t)*(Y==1))/n_pos
  }
  height = (beta[-1]+beta[-length(beta)])/2
  width = abs(diff(alpha)) # = diff(rev(omspec))
  area = sum(height*width)
  return(list(area=area,alpha=alpha,beta=beta))
}

bininv <- function(m, e, c){
  plo = 0
  phi = 1
  p = .5
  
  max_iter = 20
  tol = c*.001
  
  iter = 0
  while(iter <= max_iter){
    iter = iter + 1
    bintail = pbinom(e, m, p)
    if(abs(bintail-c)<=tol){
      return(p)
    }
    if(bintail<c){
      phi = p
    }else{
      plo = p
    }
    p = (phi+plo)/2
  }
  return(p)
}

kappahat <- function(posteriors, Y, n1, n2){
  ## Extract the empirical false and true positive probabilities
  nauc = 200 ## number of sections when calculating AUC
  auc.fit = auc(posteriors, Y, nauc)
  a = auc.fit$area
  fp = auc.fit$alpha
  tp = auc.fit$beta
  ind = (fp < 1) & (fp > 0) ## omit endpoints if present
  false_positive_rates = fp[ind]
  true_positive_rates = tp[ind]
  n_roc = length(false_positive_rates)
  
  ## Estimate slopes using binomial tail inversion
  numer = rep(0, n_roc)
  denom = rep(0, n_roc)
  delta = .1 ## confidence level
  for(ii in 1:n_roc){
    numer[ii] = bininv(n2, n2*(1 - true_positive_rates[ii]), delta)
    denom[ii] = 1 - bininv(n1, n1*false_positive_rates[ii], delta)
  }
  
  ## take minimum estimated slope
  slopes = numer / denom
  s = min(slopes)
  
  return(s)
}

cpe <- function(sample1, sample0, train_frac=.2, seed = 2017){
  set.seed(seed)
  sample1 = as.matrix(sample1)
  sample0 = as.matrix(sample0)
  p = ncol(sample1)
  n1 = nrow(sample1)
  n0 = nrow(sample0)
  n = n1 + n0
  X = rbind(sample0, sample1)
  Y = c(-rep(1,n0), rep(1,n1)) 
  
  perm = sample(n)
  X = X[perm,]
  Y = Y[perm]
  
  ## set up data for hold out estimation
  a_idx = 1:floor(n*train_frac)
  b_idx = setdiff(1:n, a_idx) ## estimate class probabilities on this set
  Xa = X[a_idx,]
  Xb = X[b_idx,]
  na = length(a_idx)
  nb = length(b_idx)
  Ya = Y[a_idx]
  Yb = Y[b_idx]
  n0a = sum(Ya == -1)
  n1a = sum(Ya == 1)
  n0b = sum(Yb == -1)
  n1b = sum(Yb == 1)
  nvec = c(n0a,n1a,n0b,n1b)
  
  ## jaakkola heuristic for Gaussian kernel bandwidth
  dist_mat = as.matrix(pdist(Xa[Ya==-1,], Xa[Ya==1,]))
  jaakk_heur = median(dist_mat)
  
  ## sample random Fourier directions and angles
  D = 500 ## RFF dimension
  Omega = matrix(rnorm(p*D), nrow = p) ## RVs defining RFF tranform
  beta = runif(D)*2*pi
  nauc = 20 ## number of sections when calculating AUC
  
  ## select parameters using cross-validation, maximizing AUC
  best_auc = 0
  for(sigma in jaakk_heur*2^seq(-3,1,length.out = 5)){ ## kernel bandwidth
    Z = transformFeatures(Xa/sigma, Omega, beta) ## calculate RFFs
    for(lambda in 10^seq(-2,2,length.out = 5)){
      cv_idx = cut(sample(na), breaks = 5, labels = F)
      cv_auc = rep(0,5)
      for(ii in 1:3){
        Y_train = Ya[cv_idx != ii]
        Z_train = Z[cv_idx != ii,]
        Y_test = Ya[cv_idx == ii]
        Z_test = Z[cv_idx == ii,]
        model = LiblineaR(Z_train, Y_train, type = 7, cost = 1/lambda)
        fit = predict(model, Z_test, proba = T)
        eta = fit$probabilities[,1]
        auc.fit = auc(eta, Y_test)
        cv_auc[ii] = auc.fit$area
      }
      mean_auc = mean(cv_auc)
      if(!is.na(mean_auc) & mean_auc > best_auc){
        best_auc = mean_auc
        best_params = c(lambda, sigma)
      }
    }
  }
  
  ## train on the whole set
  lambda = best_params[1]
  sigma = best_params[2]
  Za = transformFeatures(Xa/sigma, Omega, beta)
  model = LiblineaR(Za, Ya, type = 7, cost = 1/lambda)
  
  Z = transformFeatures(X/sigma, Omega, beta)
  Y_fit =  predict(model, Z, proba = T)
  
  ## calculate prediction for held out data
  Zb = transformFeatures(Xb/sigma, Omega, beta)
  fit = predict(model, Zb, proba = T)
  eta = fit$probabilities[,1]
  labels = Yb
  
  return(list(eta=eta,labels=labels,nvec=nvec,best_params=best_params, Y_fit=Y_fit))
}

# eta = c(.2, .3, .7, .8, .5, .3, .4)
# Y = c(-1,1,-1,1, 1, -1, -1)
# auc(eta, Y, 3)

# X = matrix(c(1,3,1,8),ncol=2)
# Omega = matrix(c(2,3,5,8),ncol=2)
# beta = c(1,2)
# transformFeatures(X, Omega, beta)

# X = matrix(c(1,3,1,8),ncol=2)
# Omega = matrix(c(2,3,5,8),ncol=2)
# as.matrix(pdist(X, Omega))

# posteriors = c(.2, .3, .7, .8, .5, .3, .4)
# Y = c(-1,1,-1,1, 1, -1, -1)
# n1 = 3
# n2 = 4
# kappahat(posteriors, Y, n1, n2)
