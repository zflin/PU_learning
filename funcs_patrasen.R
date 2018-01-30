library(Iso)

PlotDensities <- function(a,scale=rep(1,length(a)),bw="nrd0",xlab="x"){
    for(ii in 1:length(a)){
        if(length(a[[ii]]) > 10000) a[[ii]] <- sample(a[[ii]],10000)
    }
    ds <- lapply(a,function(x){density(x,bw=bw)})
    for(ii in 1:length(ds)){
        ds[[ii]]$y <- scale[ii]*ds[[ii]]$y
    }
    xlim <- range(vapply(ds,function(x){range(x$x)},c(0,0)))
    ylim <- range(vapply(ds,function(x){range(x$y)},c(0,0)))
    par(mar=c(5,5,1,1))
    plot(0,0,xlim=xlim,ylim=ylim,col=0,xlab=xlab,ylab="Density")
    for(ii in 1:length(ds)){
        points(ds[[ii]]$x,ds[[ii]]$y,type='l',col=ii)
    }
    legend("topright",names(a),col=1:length(a),lty=1)
}


### run non-decreasing least squares
EstMixMdl <- function(data,FbD,alpha_grid=(1:200)/200){
    ## Length of the data set
    data <- sort(data)
    ## Sorts the data set
    data.1 <- unique(data)
    ## Finds the unique data points
    Fn <- ecdf(data)
    ## Computes the empirical DF of the data
    Fn.1 <- Fn(data.1)
    ## Empirical DF of the data at the data points
    ## Calculate the known F_b at the data points
    Fb <- FbD(data.1)
    ## Compute the weights (= frequency/n) of the unique data values, i.e., dF_n
    Freq <- diff(c(0,Fn.1))
    distance <- rep(0,length(alpha_grid))
    F.isS <- matrix(0,nrow=length(alpha_grid),ncol=length(Fn.1))
    for(ii in 1:length(alpha_grid)){
        ## Assumes a value of the mixing proportion
        F.hat <- (Fn.1-(1-alpha_grid[ii])*Fb)/alpha_grid[ii]
        ## Computes the naive estimator of F_s
        F.is <- pava(F.hat,Freq,decreasing=FALSE) ## Computes the Isotonic Estimator of F_s
        F.is[F.is<=0] <- 0
        F.is[F.is>=1] <- 1
        F.isS[ii,] <- F.is
        distance[ii] <- alpha_grid[ii]*sqrt(sum(((F.hat-F.is)^2)*Freq))
    }
    return(list(distance=distance,alpha_grid=alpha_grid,fs=F.isS))
}

ComputeSecondDeriv <- function(x){
    return(c(0,0,diff(diff(x))))
}

## f <- density(x)
## CreateApproxFun creates callable density
## with nice edge effects (goes to 0)
CreateApproxFun <- function(f){
    x <- f$x
    y <- f$y
    n <- length(x)
    x_min <- x[1]-(x[2]-x[1])
    x_max <- x[n] + (x[n]-x[n-1])
    x <- c(x_min,x,x_max)
    y <- c(0,y,0)
    return(approxfun(x,y,rule=2))
}

EstimateParams <- function(data,FbD,alpha_grid=(1:200)/200){
    data <- sort(data)
    out <- EstMixMdl(data,FbD,alpha_grid)
    ix <- which.max(ComputeSecondDeriv(out$distance))
    alpha_hat <- out$alpha_grid[ix]
    Fs_est <- out$fs[ix,]
    xs <- data[diff(c(0,Fs_est)) > 0]
    fs <- density(xs)
    fs <- CreateApproxFun(fs)
    f <- density(data)
    f <- CreateApproxFun(f)
    return(list(out=out,alpha_hat=alpha_hat,fs=fs,f=f,xs=xs))
}


## training data from two classes, make kde for each class
## and return classifications
OneDKDEClassifier <- function(train_cl,train_x,test_x){
    fb <- density(train_x[train_cl==1])
    fb <- CreateApproxFun(fb)
    fs <- density(train_x[train_cl==2])
    fs <- CreateApproxFun(fs)
    alpha <- mean(train_cl==2)
    num <- (1-alpha)*fb(test_x)
    den <- num + alpha*fs(test_x)
    return(num/den)
}


ComputeROCCurve <- function(cl,probs){
    cl_pred <- cl[order(probs,runif(length(probs)),decreasing=TRUE)]
    return(cbind(cumsum(cl_pred)/sum(cl==1),cumsum(-cl_pred + 1)/sum(cl==0)))
}
    

