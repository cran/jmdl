#' @title Extract or Get Generalized Components from a Fitted Joint Mean
#' Correlation Model
#'
#' @description Extract (or "get") "components" - in a generalized sense - from
#' a fitted joint mean correlation model from an object of class "JmdlMod".
#'
#' @param object a fitted joint mean correlation model of class "JmdlMod", i.e.,
#' typically the result of jmdl().
#' @param name a character vector specifying the name(s) of the "component".
#'
#' possible values are:
#' \describe{
#'   \item{\code{"m"}}{a vector of number of measurement for each subject}
#'   \item{\code{"Y"}}{response matrix}
#'   \item{\code{"X"}}{model matrix for mean structure}
#'   \item{\code{"W"}}{model matrix for correlation structure (the lower
#'   triangular matrix)}
#'   \item{\code{"offset"}}{a vecter to be added to a linear predictor}
#'   \item{\code{"theta"}}{parameter estimates of joint mean correlation model}
#'   \item{\code{"beta"}}{parameter estimates for mean structure model}
#'   \item{\code{"delta"}}{parameter estimates for mean structure model (for
#'   Nbinom model)}
#'   \item{\code{"gamma"}}{parameter estimates for correlation structure (the
#'   lower triangular matrix)}
#'   \item{\code{"stdbeta"}}{standard error for parameter beta}
#'   \item{\code{"stddelta"}}{standard error for parameter delta}
#'   \item{\code{"stdgamma"}}{standard error for parameter gamma}
#'   \item{\code{"loglik"}}{log-likelihood, except for a constant}
#'   \item{\code{"family"}}{the marginal distributions of the discrete variables}
#'   \item{\code{"q"}}{degree of polynomial of the time lag to model the lower
#'   triangular matrix}
#'   \item{\code{"time"}}{a vector of time from the data}
#' }
#'
#'
#' @examples
#'
#' mydat <- toydata
#' fit <- jmdl(Y|id|time ~ X, data = mydat, q = 2, family ='Bernoulli')
#' beta <- getJMDL(fit, "beta")
#' beta
#' loglik  <- getJMDL(fit, "loglik")
#' loglik
#'
#' @export
getJMDL <- function(object, name) UseMethod("getJMDL")


#' @describeIn getJMDL Extract or Get Generalized Components from a Fitted Joint
#' Mean Correlation Model
#' @export
getJMDL.jmdlMod <- function(object,
                            name = c("m", "Y", "X", "W", "offset", "theta", "beta", "gamma", "delta",
                                     "loglik", "family", "q", "time", "stdbeta", "stdgamma",
                                     "stddelta"))
{
  if(missing(name)) stop("'name' must not be missing")

  opt     <- object@opt
  args    <- object@args
  devcomp <- object@devcomp
  std     <- object@std
  offset  <- object@offset

  m = args$m
  Y = args$Y
  X = args$X
  W = args$W
  time = args$time
  theta  = drop(opt$par)

  if (devcomp$dims['Bernoulli']) family <- 'Bernoulli'
  if (devcomp$dims['Poisson'])   family <- 'Poisson'
  if (devcomp$dims['Nbinom'])    family <- 'Nbinom'

  switch(name,
         "m" = args$m,
         "Y" = args$Y,
         "X" = args$X,
         "W" = args$W,
         "offset" = offset,
         "time"   = args$time,
         "theta"  = drop(opt$theta),
         "beta"   = drop(opt$beta),
         "gamma"  = drop(opt$gamma),
         "delta"  = drop(opt$delta),
         "loglik" = opt$loglik,
         "q"      = object$q,
         "family" = family,
         "stdbeta"  = std$stdbeta,
         "stdgamma" = std$stdgamma,
         "stddelta" = std$stddelta)
}




##--------------------------------------------------------------------------
lagseq <- function(time)
{
  res <- NULL
  if(length(time) != 1) {
    for(i in 2:length(time)) {
      for(j in 1:(i-1))
        res <- c(res, (time[i] - time[j]))
    }
  }
  res
}




#' @title Plot Fitted Results and Model Diagnostics
#'
#' @description plot (a) fitted angles; (b) fitted correlations vs time lag;
#' (c) the empirical distribution function vs the fitted distribution function;
#' (d) the empirical correlations vs the fitted correlations, when the descrete
#' longitudinal dataset is balanced.
#'
#' @param object a fitted joint mean correlation model of class "JmdlMod", i.e.,
#' typically the result of jmdl().
#' @param time a vector of obeservation time points
#'
#'
#'
#' @export
angle.plot <- function(object, time)
{
  debug <- 0
  op <- par(mfrow = c(2, 2))

  opt <- object@opt
  gamma  <- opt$gamma
  beta   <- opt$beta
  lgma <- length(gamma)

  args   <- object@args
  dims   <- object@devcomp$dims
  family = getJMDL(object,"family")

  m <- args[["m"]]
  Y <- args[["Y"]]
  X <- args[["X"]]
  W <- args[["W"]]
  W1 <- W[1,,,]
  offset <- getJMDL(object,"offset")

  if (length(unique(m)) != 1)
    stop("Unbalanced longitudinal dataset.")

  if(is.null(offset)) offset <- c(rep(0,length(m)*m[1]))

  DataMat <- Y

  tlag  <- lagseq(time)
  tslag <- seq(min(tlag), max(tlag), length.out = 100)

  W.tslag <- NULL
  for(i in 0:(lgma-1)) W.tslag <- cbind(W.tslag, tslag^i)

  Wgma <- W.tslag %*% gamma

  n <- length(m)
  zscore<-matrix(0,n,m)

  for(i in 1:n){
    if(i==1){
      xi<-X[1:m[1],,drop=F]
      osi<-offset[1:m[1]]
    }else{
      xi<-X[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i]),,drop=F]
      osi<-offset[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i])]
    }
    yi = Y[i,]

    if(family == 'Nbinom'){
      delta = opt$delta
      mui<-as.vector(exp(osi+xi%*%beta))
      sizei<-delta
      probi<-sizei/(sizei+mui)
      zscore[i,]<-qnorm(pnbinom(yi,sizei,probi))
    }
    if(family == 'Bernoulli'){
      pij<-as.vector(boot::inv.logit(osi+xi%*%beta))
      zscore[i,]<-qnorm(pbinom(yi,size=1,prob=pij))
    }
    if(family == 'Poisson'){
      phi.theta<-as.vector(exp(osi+xi%*%beta))
      zscore[i,]<-qnorm(ppois(yi,phi.theta))
    }
  }

  R <- cor(zscore)

  # first plot
  B <- t(chol(R))
  PhiMat <- matrix(0, dim(B)[1], dim(B)[2])
  for(j in 2:dim(B)[1]) {
    for(k in 1:(j-1)) {
      tmp <- 1
      if (k != 1) {
        tmp <- prod(sin(PhiMat[j, 1:(k-1)]))
      } # if
      PhiMat[j,k] <- acos(B[j, k]/tmp)
    } # for k
  } # for j
  PhiMatt <- t(PhiMat)

  phi <- PhiMatt[upper.tri(PhiMatt, diag=FALSE)]
  plot(tlag, tan(pi/2-phi), xlab="Lag", ylab=expression(tan(pi/2-w)), pch=16,main="(a)")
  lines(tslag, Wgma)

  # second plot
  fit.Tmat<-fit.phimat<-diag(1, dim(B)[1], dim(B)[2])
  for(j in 2:dim(B)[1]){
    for(k in 1:(j-1)){
      fit.phimat[j,k]<-pi/2-atan(W1[j,k,]%*%gamma)
    }
    fit.Tmat[j,j]<-prod(sin(fit.phimat[j,1:(j-1)]))
    fit.Tmat[j,1]<-cos(fit.phimat[j,1])
    if(j>2){
      for(l in 2:(j-1))
        fit.Tmat[j,l]<-cos(fit.phimat[j,l])*prod(sin(fit.phimat[j,1:(l-1)]))
    }
  }
  fit.cormat<-fit.Tmat%*%t(fit.Tmat)

  cor <- fit.cormat[upper.tri(fit.cormat, diag=FALSE)]
  plot(tlag, cor, xlab="Lag", ylab="Correlation",pch=16,main="(b)")

  # third plot
  plot(R[lower.tri(cor(zscore))],
       fit.cormat[lower.tri(fit.cormat)],
       xlim=c(0.2,0.8),ylim=c(0.2,0.8),
       xlab="Sample-based Correlation",
       ylab="Fitted Correlation",pch=16,main="(c)")
  abline(0,1)


  # forth plot
  Ecdf = c(0)
  for(i in 1:n){
    sum = 0
    for(j in 1:n){
      temp = prod(I(Y[j,]<Y[i,]))
      sum  = sum + temp
    }
    Ecdf[i] = sum/n
  }
  f.fit <- c(0)
  for(i in 1:n){
    f.fit[i] <- mvtnorm::pmvnorm(c(rep(-Inf, m[i])),zscore[i,],corr=fit.cormat)
  }
  plot(sort(f.fit),sort(Ecdf),
       xlab="Fitted CDF",
       ylab="Empirical CDF",pch=16,main="(d)")
  abline(0,1)
}




#' @title Plot Fitted Results
#'
#' @description plot (a) fitted angles; (b) fitted correlations vs time lag, when
#' the descrete longitudinal dataset is unbalanced.
#'
#' @param object a fitted joint mean correlation model of class "JmdlMod", i.e.,
#' typically the result of jmdl().
#'
#'
#' @export
un.angle.plot <- function(object)
{
  debug <- 0
  op <- par(mfrow = c(1, 2))

  opt <- object@opt
  gamma  <- opt$gamma
  beta   <- opt$beta
  lgma <- length(gamma)

  args   <- object@args
  dims   <- object@devcomp$dims
  family = getJMDL(object,"family")
  stdgma = getJMDL(object,"stdgamma")

  m <- args[["m"]]
  Y <- args[["Y"]]
  X <- args[["X"]]
  W <- args[["W"]]


  DataMat <- Y

  time <- args[["time"]]

  n <- length(m)
  ncor <- NULL
  ntlag <- NULL


  # first plot
  for(i in 1:n){
    tlag = NULL
    cor = NULL
    fit.Tmat<-fit.phimat<-diag(1, m[i], m[i])
    if(m[i]>1){
      for(j in 2:m[i]){
        for(k in 1:(j-1)){
          fit.phimat[j,k]<-pi/2-atan(W[i,j,k,]%*%gamma)
        }
        fit.Tmat[j,j]<-prod(sin(fit.phimat[j,1:(j-1)]))
        fit.Tmat[j,1]<-cos(fit.phimat[j,1])
        if(j>2){
          for(l in 2:(j-1))
            fit.Tmat[j,l]<-cos(fit.phimat[j,l])*prod(sin(fit.phimat[j,1:(l-1)]))
        }
      }
      if(i==1){
        timei<-time[1:m[1],drop=F]
      }else{
        timei<-time[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i]),drop=F]
      }
      tlag <- lagseq(timei)
    }
    fit.cormat <- fit.Tmat%*%t(fit.Tmat)
    cor <- fit.cormat[upper.tri(fit.cormat, diag=FALSE)]
    ncor <- c(ncor, cor)
    ntlag <- c(ntlag, tlag)
  }
  plot(ntlag, ncor, xlab="Lag", ylab="Correlation",pch=20,cex=0.5,main="(a)")


  # second plot
  tslag <- seq(0, max(time) - min(time), length.out = n)
  W.tslag <- NULL

  for(i in 0:(lgma-1)) W.tslag <- cbind(W.tslag, tslag^i)
  Wgma <- drop(W.tslag %*% gamma)

  ft.sd<-(W.tslag%*%stdgma)

  plot(tslag, Wgma, type = 'l', xlab = "Lag", ylab = expression(tan(pi/2-w)), ylim=c(-1,2),main="(b)")
  lines(tslag, Wgma+1.96*ft.sd, lty = 2, lwd = 2)
  lines(tslag, Wgma-1.96*ft.sd, lty = 2, lwd = 2)
}



#' @title Pairwise Likelihood Ratio Statistic Test
#'
#' @description Conducts a pairwise likelihood ratio test for joint mean-correlation regression.
#'
#' @param fit a fitted joint mean correlation model of class "JmdlMod", i.e.,
#' typically the result of jmdl().
#' @param id the id of paremeter to test
#'
#'
#' @export
lrt.test <- function(fit, id){
  #theta0 = 0
  m <- getJMDL(fit, "m")
  Y <- getJMDL(fit, "Y")
  X <- getJMDL(fit, "X")
  W <- getJMDL(fit, "W")
  family <- getJMDL(fit, "family")
  offset <- getJMDL(fit, "offset")
  lrt.full <- fit@opt
  id <- c(id)
  theta0 <- c(rep(0,length(id)))

  lgma<-dim(W)[4]
  Yv<-na.exclude(as.vector(t(Y)))

  ##Bernoulli
  if(family == 'Bernoulli'){
    #null comlik estimation
    y.fit0<-glm(Yv~offset(offset)+X[,-id]-1,family='binomial')
    theta<-c(coef(y.fit0),rep(0,lgma))
    lrt0<-minqa::bobyqa(par=theta,fn=cond.comlik, theta0=theta0,id=id,Y=Y,X=X,W=W,m=m,
                        family=family, offset=offset)

    lm.obj <- glm(Yv ~ offset(offset) + X-1 ,family='binomial')
    bta1 <- coef(lm.obj)
    gma1 <- rep(0, lgma)
    theta1 <- c(bta1, gma1)
    lrt1<-minqa::bobyqa(par=theta1,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
  }

  ##Poisson
  if(family == 'Poisson'){
    #null comlik estimation
    y.fit0<-glm(Yv~offset(offset)+X[,-id]-1,family='poisson')
    theta<-c(coef(y.fit0),rep(0,lgma))
    lrt0<-minqa::bobyqa(par=theta,fn=cond.comlik, theta0=theta0, id=id,Y=Y,X=X,W=W,m=m,
                        family=family, offset=offset)

    lm.obj <- glm(Yv ~ offset(offset) + X-1 ,family='poisson')
    bta1 <- coef(lm.obj)
    gma1 <- rep(0, lgma)
    theta1 <- c(bta1, gma1)
    lrt1<-minqa::bobyqa(par=theta1,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
  }

  ##Nbinom
  if(family == 'Nbinom'){
    lbta = ncol(X)
    if(sum(id==(lbta+1))) return(list(lrt=NA, pval=NA, critval=NA))
    #null comlik estimation
    y.fit0<-MASS::glm.nb(Yv~offset(offset)+X[,-id]-1)
    theta<-c(coef(y.fit0),y.fit0$theta,rep(0,lgma))
    lrt0<-minqa::bobyqa(par=theta,fn=cond.comlik, theta0=theta0, id=id,Y=Y,X=X,W=W,m=m,
                        family=family,offset=offset)

    lm.obj <- MASS::glm.nb(Yv ~ offset(offset) + X-1)
    bta1 <- coef(lm.obj)
    gma1 <- rep(0, lgma)
    theta1 <- c(bta1, lm.obj$theta, gma1)
    lrt1<-minqa::bobyqa(par=theta1,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
  }

  lrt<-2*(lrt0$fval-lrt1$fval)

  theta<-rep(0,length(lrt1$par))
  theta[id]<-theta0
  theta[-id]<-lrt1$par[-id]

  #pvalue
  #set.seed(1234)
  cov.fit0<-asycov(m=m,Y=Y,X=X,W=W,family=family,theta=theta,offset=offset)
  chimat<-matrix(rchisq(length(id)*10000, df=1),ncol=length(id))
  if(length(id)==1){
    lmds<-cov.fit0$Gmat.inv[id,id]/cov.fit0$Kmat.inv[id,id]
    a<-lmds*chimat
  }else{
    lmds<-eigen(solve(cov.fit0$Kmat.inv[id,id])%*%cov.fit0$Gmat.inv[id,id])
    a<-chimat%*%lmds$values
  }
  a<-as.vector(a)
  pval<-sum(a>lrt)/10000
  return(list(lrt.statistic=lrt, pval=pval, critval=round(quantile(a,c(0.025, 0.975)), 4))) #up 95% quantile
}
