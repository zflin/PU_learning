#################################################################################
## this code would create two temporary data sets 
##      * waveform.RData
##      * tcdb_swissprot.RData
#################################################################################

if(!dir.exists("figs")){
  dir.create("figs")
}
if(!dir.exists("data")){
  dir.create("data")
}

#============================== data: waveform =============================
##' create waveform.RData
##' dat1 : from class == 1
##' dat2 : from class == 2

rm(list = ls())
library(mlbench)
set.seed(17)
n = 6e4
p<-mlbench.waveform(n)

ind = (p$classes == "1")
dat1 <- as.data.frame(cbind(p$x[ind,], p$classes[ind]))
names(dat1)[22] <- "cl"
dat1$cl <- as.character(dat1$cl)

ind = (p$classes == "2")
dat2 <- as.data.frame(cbind(p$x[ind,], p$classes[ind]))
names(dat2)[22] <- "cl"
dat2$cl <- as.character(dat2$cl)

save(dat1, dat2, file = 'data/waveform.RData')

#============================== data: TCDB-SwissProt ===============================

rm(list=ls())

library(devtools)

if(!dir.exists("data/tcdb_swissprot")){
  dir.create("data/tcdb_swissprot")
}
# download.file("https://raw.github.com/zflin/data/master/tcdb_swissprot/N.zip",
#               destfile = "data/tcdb_swissprot/N.zip", method = "wget")
# download.file("https://raw.github.com/zflin/data/master/tcdb_swissprot/P.zip",
#               destfile = "data/tcdb_swissprot/P.zip", method = "wget")
# download.file("https://raw.github.com/zflin/data/master/tcdb_swissprot/Q.zip",
#               destfile = "data/tcdb_swissprot/Q.zip", method = "wget")

download.file("https://github.com/zflin/PU_learning/blob/master/data/tcdb_swissprot/N.zip?raw=true",
              destfile = "data/tcdb_swissprot/N.zip", method = "wget")
download.file("https://github.com/zflin/PU_learning/blob/master/data/tcdb_swissprot/P.zip?raw=true",
              destfile = "data/tcdb_swissprot/P.zip", method = "wget")
download.file("https://github.com/zflin/PU_learning/blob/master/data/tcdb_swissprot/Q.zip?raw=true",
              destfile = "data/tcdb_swissprot/Q.zip", method = "wget")


unzip("data/tcdb_swissprot/N.zip", exdir = "data/tcdb_swissprot/")
unzip("data/tcdb_swissprot/P.zip", exdir = "data/tcdb_swissprot/")
unzip("data/tcdb_swissprot/Q.zip", exdir = "data/tcdb_swissprot/")

fs_P = list.files("data/tcdb_swissprot/P")
fs_N = list.files("data/tcdb_swissprot/N")
fs_Q = list.files("data/tcdb_swissprot/Q")

P = list()
for(ii in 1:length(fs_P)){
  #ii = 33
  x = readLines(paste0('data/tcdb_swissprot/P/',fs_P[ii]))
  tmp = cumsum(unname(unlist(sapply(x, function(y){startsWith(y," ") | startsWith(y,"-!-") | startsWith(y,"\"") | endsWith(y,"\";")}))))
  tmp1 = x[which(tmp==max(tmp))[-1]]
  tmp2 = paste(tmp1,collapse = " ")
  tmp3 = unlist(strsplit(tmp2, split = ". ", fixed = T))
  tmp4 = c(unlist(strsplit(tmp3[1], split = "; ")),tmp3[2])
  P[[ii]] = tmp4
}
N = list()
for(ii in 1:length(fs_N)){
  #ii = 33
  x = readLines(paste0('data/tcdb_swissprot/N/',fs_N[ii]))
  tmp = cumsum(unname(unlist(sapply(x, function(y){startsWith(y," ") | startsWith(y,"-!-") | startsWith(y,"\"") | endsWith(y,"\";")}))))
  tmp1 = x[which(tmp==max(tmp))[-1]]
  tmp2 = paste(tmp1,collapse = " ")
  tmp3 = unlist(strsplit(tmp2, split = ". ", fixed = T))
  tmp4 = c(unlist(strsplit(tmp3[1], split = "; ")),tmp3[2])
  N[[ii]] = tmp4
}
Q = list()
for(ii in 1:length(fs_Q)){
  #ii = 33
  x = readLines(paste0('data/tcdb_swissprot/Q/',fs_Q[ii]))
  tmp = cumsum(unname(unlist(sapply(x, function(y){startsWith(y," ") | startsWith(y,"-!-") | startsWith(y,"\"") | endsWith(y,"\";")}))))
  tmp1 = x[which(tmp==max(tmp))[-1]]
  tmp2 = paste(tmp1,collapse = " ")
  tmp3 = unlist(strsplit(tmp2, split = ". ", fixed = T))
  tmp4 = c(unlist(strsplit(tmp3[1], split = "; ")),tmp3[2])
  Q[[ii]] = tmp4
}


tmp1 = unlist(P)
tmp2 = unlist(N)
tmp3 = unlist(Q)
feat = unique(c(tmp1,tmp2,tmp3))

n1 = length(P)
n2 = length(N)
n3 = length(Q)
n = n1 + n2 + n3
p = length(feat)
dat0 = matrix(0,nrow = n, ncol = p)
for(ii in 1:n1){
  dat0[ii,feat %in% P[[ii]]] = 1
}
for(ii in 1:n2){
  dat0[ii+n1,feat %in% N[[ii]]] = 1
}
for(ii in 1:n3){
  dat0[ii+n1+n2,feat %in% Q[[ii]]] = 1
}

dat0_pca = princomp(dat0)
cum_rate = cumsum(dat0_pca$sdev^2)/sum(dat0_pca$sdev^2)
#plot(tmp1)
k = 200
cum_rate[k]
dat = dat0_pca$scores[,1:k]
dat = as.data.frame(dat)
dat$cl = as.factor(c(rep("labeled",n1),rep("unlabeled",n2+n3)))

alpha = n2/(n2+n3)
cl = c(rep(0,n2),rep(1,n3))

save(dat, alpha, cl, file = "data/tcdb_swissprot.RData")

