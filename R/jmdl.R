#' @import stats
NULL

#' @import graphics
NULL

#' epilepsy2 Data
#'
#' Data from a clinical trial of 59 epileptics. For a baseline, patients were observed for 8 weeks
#' and the number of seizures recorded. The patients were then randomized to treatment by the drug
#' Progabide (31 patients) or to the placebo group (28 patients). They were observed for four
#' 2-week periods and the number of seizures recorded.
#'
#' \itemize{
#'   \item id: identifying number
#'   \item treat: 1=treated, 0=not
#'   \item seizures: number of seizures
#'   \item expind: 0=baseline period, 1=treatment period
#'   \item timeadj: weeks of period
#'   \item age: in years
#'   \item time: in weeks
#' }
#'
#' @docType data
#' @keywords datasets
#' @name epilepsy2
#' @usage data(epilepsy2)
#' @format A data frame with 280 rows and 7 variables
NULL


#' pbc3 Data
#'
#' Followup of 312 randomised patients with primary biliary cirrhosis, a rare autoimmune
#' liver disease, at Mayo Clinic.
#'
#'
#' \itemize{
#'   \item id: patients identifier; in total there are 312 patients
#'   \item years: number of years between registration and the earlier of death, transplantion, or study analysis time
#'   \item status: a factor with levels alive, transplanted and dead
#'   \item drug: a factor with levels placebo and D-penicil
#'   \item age: at registration in years
#'   \item sex: a factor with levels male and female
#'   \item year: number of years between enrollment and this visit date, remaining values on the line of data refer to this visit
#'   \item ascites: a factor with levels No and Yes
#'   \item hepatomegaly: a factor with levels No and Yes
#'   \item spiders: a factor with levels No and Yes
#'   \item edema: a factor with levels No edema (i.e., no edema and no diuretic therapy for edema), edema no diuretics (i.e., edema present without diuretics, or edema resolved by diuretics), and edema despite diuretics (i.e., edema despite diuretic therapy)
#'   \item serBilir: serum bilirubin in mg/dl
#'   \item albumin: albumin in gm/dl
#'   \item platelets: platelets per cubic ml / 1000
#'   \item prothrombin: prothrombin time in seconds
#' }

#'
#' @docType data
#' @keywords datasets
#' @name pbc3
#' @usage data(pbc3)
#' @format A data frame with 1390 rows and 15 variables
NULL


#' toy Data
#'
#' A simulation dataset.
#'
#' \itemize{
#'   \item id: identifying number
#'   \item Y: response variable: 1, 0
#'   \item time: time
#'   \item X: independent variable
#' }
#'
#' @docType data
#' @keywords datasets
#' @name toydata
#' @usage data(toydata)
#' @format A data frame with 50 rows and 4 variables
NULL


#' @title Fit Joint Mean-Correlation Models For Discrete Longitudinal Data
#'
#' @description Fit a joint mean-correlation model to discrete longitudinal data.
#'
#' @param formula a two-sided linear formula object describing the correlation
#' for both the mean and correlation matrix part of the model, with the response,
#' the corresponding subject id and measurement time on the left of a operator~,
#' divided by vertical bars ("|").
#' @param data data frame containing the variables named in formula.
#' @param q degree of polynomial of the time lag to model the lower triangular matrix.
#' @param W.appendix appendix array to model time-dependent covariates for the lower triangular matrix.
#' @param offset a term to be added to the linear predictor.
#' @param theta starting values for the parameters in the model.
#' @param family  the marginal distributions of the discrete variables.
#' choose 'Bernoulli', 'Poisson' or 'Nbinom'.
#'
#' @examples
#' data(toydata)
#' mydat <- toydata
#' fit <- jmdl(Y|id|time ~ X, data = mydat, q = 2, family ='Bernoulli')
#' @export
jmdl <- function(formula, data = NULL, q = 2, theta = NULL, W.appendix = NULL,
                 offset = NULL, family = c('Bernoulli', 'Nbinom', 'Poisson'))
{
  mc <- mcout <- match.call()
  if(missing(family))
    stop("method must be specified")

  mc[[1]] <- quote(ldFormula)
  args <- eval(mc, parent.frame(1L))

  opt <- do.call(OptimizeJmdl,
                 c(args, list(family=family),list(offset=offset),list(theta=theta)))

  lth <- length(opt$par)
  lb  <- dim(args$X)[2]
  if(is.null(W.appendix)){
    qq = q
  }else{
    qq = q + dim(W.appendix)[4]
  }

  lopt<-list(beta=opt$par[1:lb] , delta=opt$par[lb+1] ,
             gamma=opt$par[(lth-qq):lth], theta = opt$par,
             loglik=opt$fval )

  if(is.null(offset)) offset = c(rep(0,length(args$time)))

  asy <- asycov(args$m, args$Y, args$X, args$W, family=family, theta = opt$par, offset = offset)
  std <- sqrt(diag(asy$Gmat.inv))
  lstd <- list(stdbeta=std[1:lb] , stddelta=std[lb+1] ,
               stdgamma=std[(lth-qq):lth])

  tval  <- NULL
  ltval <- list(tval)
  p    <- NULL
  lp   <- list(p)

  for(i in 1:lth){
    ptest <- p.test(c(i), c(0), opt, std, family=family)
    tval[i] <- ptest$tval
    p[i]  <- ptest$pval
  }
  ltval <- list(tbeta=tval[1:lb], tdelta=tval[lb+1], tgamma=tval[(lth-qq):lth])
  lp <- list(pbeta=p[1:lb], pdelta=p[lb+1], pgamma=p[(lth-qq):lth])


  JmdlMod(opt=lopt, args=args, std = lstd, tval=ltval, p=lp, q=q, family=family, offset=offset, mc=mcout)
}

#' @title Modular Functions for Joint Mean Correlation Model Fits
#'
#' @description Modular Functions for joint mean correlation model fits
#'
#' @param formula a two-sided linear formula object describing the covariates
#' for both the mean and correlation matrix part of the model, with the response,
#' the corresponding subject id and measurement time on the left of a operator~,
#' divided by vertical bars ("|").
#' @param data a data frame containing the variables named in formula.
#' @param q degree of polynomial of the time lag to model the lower triangular matrix.
#' @param W.appendix appendix array to model time-dependent covariates for the lower triangular matrix.
#' @param theta starting values for the parameters in the model.
#' @param offset a term to be added to a linear predictor.
#' @param family  the marginal distributions of the discrete variables.
#' choose 'Bernoulli', 'Poisson' or 'Nbinom'.
#' @param m an integer vector of number of measurements for each subject.
#' @param Y a matrix of responses for all subjects.
#' @param X model matrix for mean structure model.
#' @param W model array for the lower triangular matrix.
#' @param time a vector of time from the data.
#' @param opt optimized results returned by OptimizJmdl.
#' @param args arguments returned by ldFormula.
#' @param std standard error for parameters.
#' @param tval t statistic.
#' @param p p.value.
#' @param mc matched call from the calling function.
#'
#' @name modular
NULL
#> NULL

#' @rdname modular
#' @export
ldFormula <- function(formula, data = NULL, q = 2, theta = NULL, W.appendix = NULL,
                      offset = NULL, family = c('Bernoulli', 'Nbinom', 'Poisson'))
{
  mf <- mc <- match.call()
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1, m)]

  f <- Formula::Formula(formula)
  mf[[1]] <- as.name("model.frame")
  mf$formula <- f
  mf <- eval(mf, parent.frame())

  Y    <- Formula::model.part(f, data = mf, lhs = 1)
  id   <- Formula::model.part(f, data = mf, lhs = 2)
  time <- Formula::model.part(f, data = mf, lhs = 3)

  X <- model.matrix(f, data = mf)

  index <- order(id, time)

  Y0    <- Y[index, ]
  id   <- id[index, ]
  time <- time[index, ]

  m <- table(id)
  attr(m, "dimnames") <- NULL
  n <- length(m)

  X <- X[index, ]

  if(is.null(W.appendix)){
    W<-array(NA,dim=c(n,max(m),max(m)-1,q+1))
    for(i in 1:n){
      if(i==1) {
        ti<-time[1:m[1]]
      }else{
        ti<-time[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i])]
      }

      if(m[i]>1){
        for(j in 2:m[i]){
          for(k in 1:(j-1))
            W[i,j,k,]<-c((ti[j]-ti[k])^(0:q), W.appendix[i,j,k,])
        }
      }
    }
  }else{
    h=dim(W.appendix)[4]
    W<-array(NA,dim=c(n,max(m),max(m)-1,q+1+h))
    for(i in 1:n){
      if(i==1) {
        ti<-time[1:m[1]]
      }else{
        ti<-time[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i])]
      }

      if(m[i]>1){
        for(j in 2:m[i]){
          for(k in 1:(j-1))
            W[i,j,k,]<-c((ti[j]-ti[k])^(0:q), W.appendix[i,j,k,])
        }
      }
    }
  }

  Y<-array(NA,dim=c(n,max(m)))
  for(i in 1:n){
    if(i==1){
      Y[i,1:m[i]]<-Y0[1:m[1]]
    }else{
      Y[i,1:m[i]]<-Y0[(sum(m[1:(i-1)])+1):sum(m[1:i])]
    }
  }

  list(m = m, Y = Y, X = X, W = W, time = time)
}


#' @rdname modular
#' @export
##OptimizeJmdl
OptimizeJmdl <- function(m, Y, X, W, time, offset = NULL, theta = NULL, family)
{
  missStart <- is.null(theta)

  lbta <- ncol(X)
  lgma <- dim(W)[4]

  if(is.null(offset)){
    offset = c(rep(0,length(time)))
    if (family == 'Bernoulli') {
      if (missStart) {
        Yv<-na.exclude(as.vector(t(Y)))
        lm.obj <- glm(Yv ~ X-1 ,family=binomial)
        bta0 <- coef(lm.obj)
        gma0 <- rep(0, lgma)
        theta <- c(bta0, gma0)
      }
      #est<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family)
      rlt<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      gma<-rlt$par[-c(1:lbta)]
      gma<-optim(par=gma,fn=gmalik,method="BFGS",bta.delta=rlt$par[1:lbta],m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      theta<-c(rlt$par[1:lbta],gma$par)
      est<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
    }

    if (family == 'Poisson') {
      if (missStart) {
        Yv<-na.exclude(as.vector(t(Y)))
        lm.obj <- glm(Yv ~ X-1,family='poisson')
        bta0 <- coef(lm.obj)
        gma0 <- rep(0, lgma)
        theta <- c(bta0, gma0)
      }
      rlt<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      gma<-rlt$par[-c(1:lbta)]
      gma<-optim(par=gma,fn=gmalik,method="BFGS",bta.delta=rlt$par[1:lbta],m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      theta<-c(rlt$par[1:lbta],gma$par)
      est<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
    }

    if (family == 'Nbinom') {
      if (missStart) {
        Yv<-na.exclude(as.vector(t(Y)))
        lm.obj <- MASS::glm.nb(Yv ~ X-1)
        bta0 <- coef(lm.obj)
        gma0 <- rep(0, lgma)
        theta <- c(bta0, lm.obj$theta, gma0)
      }
      #est<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family)
      rlt<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      gma<-rlt$par[-c(1:(lbta+1))]
      gma<-optim(par=gma,fn=gmalik,method="BFGS",bta.delta=rlt$par[1:(lbta+1)],m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      theta<-c(rlt$par[1:(lbta+1)],gma$par)
      est<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
    }
  }else{
    if (family == 'Bernoulli') {
      if (missStart) {
        Yv<-na.exclude(as.vector(t(Y)))
        lm.obj <- glm(Yv ~ offset(offset) + X-1 ,family=binomial)
        bta0 <- coef(lm.obj)
        gma0 <- rep(0, lgma)
        theta <- c(bta0, gma0)
      }
      #est<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family)
      rlt<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      gma<-rlt$par[-c(1:lbta)]
      gma<-optim(par=gma,fn=gmalik,method="BFGS",bta.delta=rlt$par[1:lbta],m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      theta<-c(rlt$par[1:lbta],gma$par)
      est<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
    }

    if (family == 'Poisson') {
      if (missStart) {
        Yv<-na.exclude(as.vector(t(Y)))
        lm.obj <- glm(Yv ~ offset(offset) + X-1,family='poisson')
        bta0 <- coef(lm.obj)
        gma0 <- rep(0, lgma)
        theta <- c(bta0, gma0)
      }
      rlt<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      gma<-rlt$par[-c(1:lbta)]
      gma<-optim(par=gma,fn=gmalik,method="BFGS",bta.delta=rlt$par[1:lbta],m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      theta<-c(rlt$par[1:lbta],gma$par)
      est<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
    }

    if (family == 'Nbinom') {
      if (missStart) {
        Yv<-na.exclude(as.vector(t(Y)))
        lm.obj <- MASS::glm.nb(Yv ~ offset(offset) + X-1)
        bta0 <- coef(lm.obj)
        gma0 <- rep(0, lgma)
        theta <- c(bta0, lm.obj$theta, gma0)
      }
      #est<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family)
      rlt<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      gma<-rlt$par[-c(1:(lbta+1))]
      gma<-optim(par=gma,fn=gmalik,method="BFGS",bta.delta=rlt$par[1:(lbta+1)],m=m,Y=Y,X=X,W=W,family=family,offset=offset)
      theta<-c(rlt$par[1:(lbta+1)],gma$par)
      est<-minqa::bobyqa(par=theta,fn=comlik,m=m,Y=Y,X=X,W=W,family=family,offset=offset)
    }
  }
  est
}




#' @rdname modular
#' @export
JmdlMod <- function(opt, args, std, tval, p, q, family, offset, mc)
{
  if(missing(mc))
    mc <- match.call()

  isBernoulli <- (family == "Bernoulli")
  isPoisson <- (family == "Poisson")
  isNbinom <- (family == "Nbinom")

  dims  <- c(nsub = length(args$m), max.nobs = max(args$m), q = q,
             Bernoulli = isBernoulli, Poisson = isPoisson, Nbinom = isNbinom)

  new("jmdlMod", call=mc, opt=opt, args= args, std = std, tval = tval, p = p, q=q,
      offset = offset, devcomp=list(dims=dims))
}


###----- Printing etc ----------------------------
methodTitle <- function(object, dims = object@devcomp$dims)
{
  Bernoulli <- dims[["Bernoulli"]]
  Poisson <- dims[["Poisson"]]
  Nbinom <- dims[["Nbinom"]]
  kind <- switch(Bernoulli * 1L + Poisson * 2L + Nbinom * 3L, "Bernoulli", "Poisson", "Nbinom")
  paste("Mean-Correlation Regression based on", kind)
}



cat.f <- function(...) cat(..., fill = TRUE)

.prt.methTit <- function(mtit, class) {
  cat(sprintf("%s ['%s']\n", mtit, class))
}

.prt.call <- function(call, long = TRUE) {
  if (!is.null(cc <- call$formula))
    cat.f("Formula:", deparse(cc))    ##deparse:Turn unevaluated expressions into character strings.
  if (!is.null(cc <- call$q))
    cat.f(" q:", deparse(cc))
}

.prt.loglik <- function(n2ll, digits=4)
{
  t.4 <- round(n2ll, digits)
  cat.f("logLik:", t.4)
}


print.jmdlMod <- function(x, digits = 4)
{
  dims <- x@devcomp$dims
  .prt.methTit(methodTitle(x, dims = dims), class(x))
  .prt.call(x@call); cat("\n")
  .prt.loglik(x@opt$loglik); cat("\n")


  mp = data.frame(Estimate = drop(x@opt$beta), Std.Error = drop(x@std$stdbeta),
                    t.value = drop(x@tval$tbeta), p.value = drop(x@p$pbeta) )

  row.names(mp)<- c(colnames(x@args$X))
  cat("Mean Parameters:\n")
  print(mp , digits = digits, print.gap = 2L, quote = FALSE)

  if(dims["Nbinom"]==1){
    odp = data.frame(Estimate = drop(x@opt$delta), Std.Error = drop(x@std$stddelta),
                       t.value = drop(x@tval$tdelta), p.value = drop(x@p$pdelta) )
    row.names(odp)<-c("delta")
    cat("\n Over-dispersion Parameter:\n")
    print(odp, digits = digits, print.gap = 2L, quote = FALSE)
  }

  cat("\n Correlation Parameters:\n")

  cp = data.frame(Estimate = drop(x@opt$gamma), Std.Error = drop(x@std$stdgamma),
                    t.value = drop(x@tval$tgamma), p.value = drop(x@p$pgamma) )

  row.names(cp)<- c(paste("gamma",1:length(drop(x@opt$gamma))))
  print(cp, digits = digits, print.gap = 2L, quote = FALSE)

  invisible(x)
}



#' Print information for jmdlMod-class
#'
#' @param object a fitted joint mean correlation model of class "jmdlMod", i.e.,
#' typically the result of jmdl().
#'
#' @exportMethod show
setMethod("show", "jmdlMod", function(object) print.jmdlMod(object))

summary.jmdlMod <- function(x, digits = 4)
{
  dims <- x@devcomp$dims
  .prt.methTit(methodTitle(x, dims = dims), class(x))
  .prt.call(x@call); cat("\n")
  .prt.loglik(x@opt$loglik); cat("\n")


  mp = data.frame(Estimate = drop(x@opt$beta), Std.Error = drop(x@std$stdbeta),
                  tval = drop(x@tval$tbeta), p.value = drop(x@p$pbeta) )

  row.names(mp)<- c(colnames(x@args$X))
  cat("Mean Parameters:\n")
  print(mp , digits = digits, print.gap = 2L, quote = FALSE)

  if(dims["Nbinom"]==1){
    odp = data.frame(Estimate = drop(x@opt$delta), Std.Error = drop(x@std$stddelta),
                     lrt = drop(x@tval$tdelta), p.value = drop(x@p$pdelta) )
    row.names(odp)<-c("delta")
    cat("\n Over-dispersion Parameter:\n")
    print(odp, digits = digits, print.gap = 2L, quote = FALSE)
  }

  cat("\n Correlation Parameters:\n")

  cp = data.frame(Estimate = drop(x@opt$gamma), Std.Error = drop(x@std$stdgamma),
                  lrt = drop(x@tval$tgamma), p.value = drop(x@p$pgamma) )

  row.names(cp)<- c(paste("gamma",1:length(drop(x@opt$gamma))))
  print(cp, digits = digits, print.gap = 2L, quote = FALSE)

  invisible(x)
}
