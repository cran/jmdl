
#Calculate p order polynomial
polys <- function(x,p){
  if(length(x) == 1) return(x^(0:(p-1)))
  else {
    y <- matrix(0,length(x),p)
    y[,1] <- 1
    for(i in 1:length(x))
      y[i,-1] <- x[i]^(1:(p-1))
    return(y)
  }
}


#refactor a dataframe or factor
refactor <- function (aDf){
  if(is.data.frame(aDf)){
    cat <- sapply(aDf, is.factor);
    aDf[cat] <- lapply(aDf[cat], factor);
  }else{
    aDf<-factor(aDf)
  }
  return(aDf)
}



##---------------------------------------------------------------------------
##distribution of bivnorm
dbivnorm<-function(x,y,rho){
  if(is.infinite(x)|is.infinite(y)) rlt<-0
  else rlt<-1/(2*pi*(1-rho^2))*exp(-1/(2*(1-rho^2))*(x^2-2*rho*x*y+y^2))
  return(rlt)
}


##--------------------------------------------------------------------------
##
mypnorm<-function(z1,z2,rho,sdd){
  if(!is.infinite(z1)&!is.infinite(z2))
    rlt<-pnorm((z2-z1*rho)/sdd)

  if((z2==+Inf&z1!=-Inf)|(z1==-Inf&z2!=-Inf)) rlt<-1
  if((z2==-Inf&z1!=-Inf)|(z1==+Inf&z2!=+Inf)) rlt<-0
  if(z2==-Inf & z1==-Inf) rlt<-0
  if(z2==+Inf & z1==+Inf) rlt<-0
  return(rlt)
}


##----------------------------------------------------------------------------
##derivatives of log pmf wrt to parameter bta, for m[i]=1 case
deriv.pmf <- function(bta, yi, xi, osi, family){
  if(family == "Bernoulli"){
    p <- as.vector(exp(osi+xi%*%bta)/(1+exp(osi+xi%*%bta)))
    prob <- dbinom(yi, 1, p)
    rlt <- (yi-p)*xi
  }
  if(family == "Poisson"){
    lmd <- exp(osi+xi%*%bta)
    rlt <- (yi-as.vector(lmd))*xi
  }
  return(rlt)
}



##----------------------------------------------------------------------------
##derivatives of log pmf wrt to parameter bta, for m[i]=1 case (for Nbinom)
nb.deriv.pmf<-function(bta,delta,yi,xi,osi,family){
  if(family=="Nbinom"){
    mu<- exp(sum(osi+xi*bta))
    size<-delta
    prob<-size/(size+mu)
    deriv.bta<-as.vector(yi-(size+yi)*mu/(size+mu))*xi
    deriv.delta<-digamma(yi+size)-digamma(size)-log(size+mu)-(size+yi)/(size+mu)+log(size)+1
    rlt<-c(deriv.bta,deriv.delta)
  }
  return(rlt)
}



##----------------------------------------------------------------------------
##derivative of loglikelihood function likj with respect to parameter
##(d(l_ijk)/d(theta)) (for bernoulli and poisson)
deriv.lijk<-function(j,k,bta,gma,rho,yi,xi,osi,wi,omegai,Ti,family){
  #theta=c(beta,gma)
  #rho=rhoijk
  #xi: m[i]x p covariates
  #wi: covariates for rho

  sdd<-sqrt(1-rho^2)
  Rjk<-matrix(1,2,2)
  Rjk[1,2]<-Rjk[2,1]<-rho

  if(family=="Bernoulli"){
    pij<-exp(osi[j]+xi[j,]%*%bta)/(1+exp(osi[j]+xi[j,]%*%bta))
    pik<-exp(osi[k]+xi[k,]%*%bta)/(1+exp(osi[k]+xi[k,]%*%bta))
    zij.up<-qnorm(pbinom(yi[j],1,pij))
    zij.lo<-qnorm(pbinom(yi[j]-1,1,pij))
    zik.up<-qnorm(pbinom(yi[k],1,pik))
    zik.lo<-qnorm(pbinom(yi[k]-1,1,pik))
  }

  if(family=="Poisson"){
    lmdij<-exp(osi[j]+xi[j,]%*%bta)
    lmdik<-exp(osi[k]+xi[k,]%*%bta)
    zij.up<-qnorm(ppois(yi[j],lmdij))
    zij.lo<-qnorm(ppois(yi[j]-1,lmdij))
    zik.up<-qnorm(ppois(yi[k],lmdik))
    zik.lo<-qnorm(ppois(yi[k]-1,lmdik))
  }

  #Likelihood
  Lijk<-mvtnorm::pmvnorm(c(zij.lo,zik.lo),c(zij.up,zik.up),corr=Rjk)
  if(!Lijk[1]) return(rep(0,length(bta)+length(gma)))

  #derivatives w.r.t beta
  if(zij.up==+Inf &  zik.up==+Inf)
    term1<-rep(0,length(bta))
  else
    term1<-pnorm((zik.up-rho*zij.up)/sdd)*F.deriv.bta(bta,yi[j],xi[j,],osi[j],family)+pnorm((zij.up-rho*zik.up)/sdd)*F.deriv.bta(bta,yi[k],xi[k,],osi[k],family)

  term2<-mypnorm(zij.lo,zik.up,rho,sdd)*F.deriv.bta(bta,yi[j]-1,xi[j,],osi[j],family)+mypnorm(zik.up,zij.lo,rho,sdd)*F.deriv.bta(bta,yi[k],xi[k,],osi[k],family)

  term3<-mypnorm(zij.up,zik.lo,rho,sdd)*F.deriv.bta(bta,yi[j],xi[j,],osi[j],family)+mypnorm(zik.lo,zij.up,rho,sdd)*F.deriv.bta(bta,yi[k]-1,xi[k,],osi[k],family)

  if(zik.lo==-Inf & zij.lo==-Inf)
    term4<-rep(0,length(bta))
  else
    term4<-pnorm((zik.lo-rho*zij.lo)/sdd)*F.deriv.bta(bta,yi[j]-1,xi[j,],osi[j],family)+pnorm((zij.lo-rho*zik.lo)/sdd)*F.deriv.bta(bta,yi[k]-1,xi[k,],osi[k],family)

  lijk.bta<-term1-term2-term3+term4

  #derivatives w.r.t gma
  temp<-rep(0,length(gma))
  if(j==1)
    temp<-Tijk.dev(gma,k,1,omegai,Ti,wi)
  else{
    for(s in 2:j)
      temp<-temp+ Tijk.dev(gma,j,s,omegai,Ti,wi)*Ti[k,s]+Ti[j,s]* Tijk.dev(gma,k,s,omegai,Ti,wi)
  }

  lijk.gma<-temp*(dbivnorm(zij.up,zik.up,rho)-dbivnorm(zij.lo,zik.up,rho)-dbivnorm(zij.up,zik.lo,rho)+dbivnorm(zij.lo,zik.lo,rho))
  rlt<-c(lijk.bta,lijk.gma)/Lijk

  return(rlt)
}




##------------------------------------------------------------------------------
##derivative of loglikelihood function likj with respect to parameter
##(d(l_ijk)/d(theta)) (for Nbinom)
nb.deriv.lijk<-function(j,k,bta,delta,gma,rho,yi,xi,osi,wi,omegai,Ti,family){
  #theta=c(beta,gma)
  #rho=rhoijk
  #xi: m[i]x p covariates
  #wi: covariates for rho

  sdd<-sqrt(1-rho^2)
  Rjk<-matrix(1,2,2)
  Rjk[1,2]<-Rjk[2,1]<-rho

  if(family=="Nbinom"){
    muij<-exp(osi[j]+xi[j,]%*%bta)
    muik<-exp(osi[k]+xi[k,]%*%bta)
    sizei<-delta
    probij<-sizei/(sizei+muij)
    probik<-sizei/(sizei+muik)
    zij.up<-qnorm(pnbinom(yi[j],sizei,probij))
    zij.lo<-qnorm(pnbinom(yi[j]-1,sizei,probij))
    zik.up<-qnorm(pnbinom(yi[k],sizei,probik))
    zik.lo<-qnorm(pnbinom(yi[k]-1,sizei,probik))
  }

  #Likelihood
  Lijk<-mvtnorm::pmvnorm(c(zij.lo,zik.lo),c(zij.up,zik.up),corr=Rjk)
  if(!Lijk[1]) return(rep(0,length(bta)+1+length(gma)))

  #derivatives w.r.t beta
  if(zij.up==+Inf &  zik.up==+Inf)
    term1<-rep(0,length(bta)+1)
  else
    term1<-pnorm((zik.up-rho*zij.up)/sdd)*F.deriv.mu(bta,yi[j],xi[j,],osi[j],family,delta)+pnorm((zij.up-rho*zik.up)/sdd)*F.deriv.mu(bta,yi[k],xi[k,],osi[k],family,delta)

  term2<-mypnorm(zij.lo,zik.up,rho,sdd)*F.deriv.mu(bta,yi[j]-1,xi[j,],osi[j],family,delta)+mypnorm(zik.up,zij.lo,rho,sdd)*F.deriv.mu(bta,yi[k],xi[k,],osi[k],family,delta)

  term3<-mypnorm(zij.up,zik.lo,rho,sdd)*F.deriv.mu(bta,yi[j],xi[j,],osi[j],family,delta)+mypnorm(zik.lo,zij.up,rho,sdd)*F.deriv.mu(bta,yi[k]-1,xi[k,],osi[k],family,delta)

  if(zik.lo==-Inf & zij.lo==-Inf)
    term4<-rep(0,length(bta)+1)
  else
    term4<-pnorm((zik.lo-rho*zij.lo)/sdd)*F.deriv.mu(bta,yi[j]-1,xi[j,],osi[j],family,delta)+pnorm((zij.lo-rho*zik.lo)/sdd)*F.deriv.mu(bta,yi[k]-1,xi[k,],osi[k],family,delta)

  lijk.bta<-term1-term2-term3+term4

  #derivatives w.r.t gma
  temp<-rep(0,length(gma))
  if(j==1)
    temp<-Tijk.dev(gma,k,1,omegai,Ti,wi)
  else{
    for(s in 2:j)
      temp<-temp+ Tijk.dev(gma,j,s,omegai,Ti,wi)*Ti[k,s]+Ti[j,s]* Tijk.dev(gma,k,s,omegai,Ti,wi)
  }

  lijk.gma<-temp*(dbivnorm(zij.up,zik.up,rho)-dbivnorm(zij.lo,zik.up,rho)-dbivnorm(zij.up,zik.lo,rho)+dbivnorm(zij.lo,zik.lo,rho))

  rlt<-c(lijk.bta,lijk.gma)/Lijk

  return(rlt)
}




##-----------------------------------------------------------------------------
##dF/d(beta)
##for bernoulli and poisson
F.deriv.bta<-function(bta,y,x,os,family){
  if(family=="Bernoulli"){
    #theta=logit(p)=x%*%bta+os
    theta<-x%*%bta+os
    p<-exp(theta)/(1+exp(theta))
    if(y==0)
      rlt<- -as.vector((1-p)*p)*as.vector(x)
    else
      rlt<-rep(0,length(x))
  }

  if(family=="Poisson"){
    lmd<-exp(x%*%bta+os)
    rlt<-(as.vector(ppois(y-1,lmd))-as.vector(ppois(y,lmd)))*as.vector(lmd)*as.vector(x)
  }

  return(rlt)
}



##------------------------------------------------------------------------------
##dF/d(beta)
##for Nbinom
F.deriv.mu<-function(bta,y,x,os,family,delta){
  if(family=="Nbinom"){
    mu<-exp(sum(x*bta+os))
    size<-delta
    prob<-size/(size+mu)
    rlt<-rep(0,length(bta)+1)
    for(k in 0:y){
      deriv.bta<-(k-(size+k)*mu/(size+mu))*x
      deriv.delta<-digamma(k+size)-digamma(size)-log(size+mu)-(size+k)/(size+mu)+log(size)+1
      rlt<-rlt+dnbinom(k,size,prob)*c(deriv.bta,deriv.delta)
    }
  }
  return(rlt)
}



##---------------------------------------------------------------------------
#first derivative of Tijk w.r.t gamma ( d(T_ijk/d(gamma) )
Tijk.dev<-function(gma,t,s,omegai,Ti,wi){
  #s<t
  #omega=pi/2-atan(theta)
  #theta=w%*%gma
  theta.ts<-wi[t,s,]%*%gma
  if(s>t) stop("s should be less or equal to than t")
  if(s==1)
    rlt<-as.vector(sin(omegai[t,s])* wi[t,s,])/as.vector(1+theta.ts*theta.ts)
  else{
    rlt<-as.vector(tan(omegai[t,s])*wi[t,s,])/as.vector(1+theta.ts*theta.ts)
    temp<-rep(0,dim(wi)[3])
    for(k in 1:(s-1)){
      theta.tk<-sum(wi[t,k,]*gma)
      temp<-temp-wi[t,k,]/(tan(omegai[t,k])*(1+theta.tk*theta.tk))
    }
    if(s<t)
      rlt<-Ti[t,s]*(rlt+temp)
    else
      rlt<-Ti[t,s]*temp
  }

  return(rlt)
}




##-----------------------------------------------------------------------------
##comlik
comlik <- function(theta, m, Y, X, W, family, offset)
{
  #Y: n x max(m) array
  #X: N x p covariates matrix
  #W: n x max(m) x (max(m)-1) x length(gamma) array
  #m: measurement vector of lengt n

  n<-length(m)
  p<-ncol(X)
  bta<-theta[1:p]
  gma<-theta[-c(1:p)]

  if(family == 'Nbinom'){
    bta<-theta[1:p]
    delta<-theta[p+1]
    gma<-theta[-c(1:(p+1))]
  }

  lik<-0
  for(i in 1:n){
    yi<-Y[i,1:m[i]]
    if(i==1){
      xi<-X[1:m[1],,drop=F]
      osi<-offset[1:m[1]]
    }else{
      xi<-X[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i]),,drop=F]
      osi<-offset[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i])]
    }

    if(family == 'Bernoulli'){
      pij<-as.vector(boot::inv.logit(osi + xi%*%bta))
      u.ui<-pbinom(yi,size=1,prob=pij)
      l.ui<-pbinom(yi-1,size=1,prob=pij)
    }

    if(family == 'Poisson'){
      phi.theta<-as.vector(exp(osi + xi%*%bta))
      u.ui<-ppois(yi,phi.theta)
      l.ui<-ppois(yi-1,phi.theta)
    }

    if(family == 'Nbinom'){
      mui<-as.vector(exp(osi + xi%*%bta))
      sizei<-delta
      probi<-sizei/(sizei+mui)
      u.ui<-pnbinom(yi,size=sizei,prob=probi)
      l.ui<-pnbinom(yi-1,size=sizei,prob=probi)
    }

    zi.up<-qnorm(u.ui)
    zi.lo<-qnorm(l.ui)

    Ti<-phi<-matrix(0,m[i],m[i])
    if(m[i]>1){
      Ti[1,1]<-1
      for(j in 2:m[i]){
        for(k in 1:(j-1)){
          phi[j,k]<-pi/2-atan(W[i,j,k,]%*%gma)
        }
        Ti[j,j]<-prod(sin(phi[j,1:(j-1)]))
        Ti[j,1]<-cos(phi[j,1])
        if(j>2){
          for(l in 2:(j-1))
            Ti[j,l]<-cos(phi[j,l])*prod(sin(phi[j,1:(l-1)]))
        }
      }
      Ri<-Ti%*%t(Ti)
      temp<-0
      Rkl<-diag(2)
      for(l in 2:m[i])
        for(k in 1:(l-1)){
          Rkl[1,2]<-Rkl[2,1]<-Ri[k,l]
          prob<-mnormt::biv.nt.prob(df=Inf,lower=zi.lo[c(k,l)],upper=zi.up[c(k,l)],mean=c(0,0),S=Rkl)
          if(prob<=0) return(Inf)
          temp<-temp+log(prob)
        }
    }else{
      if(family == 'Bernoulli'){
        temp<-log(dbinom(yi,size=1,prob=pij))
      }

      if(family == 'Poisson'){
        temp<-log(dpois(yi,phi.theta))
      }

      if(family == 'Nbinom'){
        temp<-log(dnbinom(yi,size=sizei,prob=probi))
      }
    }

    lik<-lik+temp
  }
  return(-lik)
}



cond.comlik <- function(theta,theta0,id, m, Y, X, W, family, offset)
{
  #Y: n x max(m) array
  #X: N x p covariates matrix
  #W: n x max(m) x (max(m)-1) x length(gamma) array
  #m: measurement vector of lengt n
  n<-length(m)
  p<-ncol(X)

  if(id <= p) {
    full.theta<-rep(0,length(theta)+length(theta0))
    full.theta[id]<-theta0
    full.theta[-id]<-theta
  }else{
    full.theta<-theta
    full.theta[id]<-theta0
  }


  bta<-full.theta[1:p]
  gma<-full.theta[-c(1:p)]

  if(family == 'Nbinom'){
    bta<-full.theta[1:p]
    delta<-full.theta[p+1]
    gma<-full.theta[-c(1:(p+1))]
  }

  lik<-0
  for(i in 1:n){
    yi<-Y[i,1:m[i]]
    if(i==1){
      xi<-X[1:m[1],,drop=F]
      osi<-offset[1:m[1]]
    }else{
      xi<-X[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i]),,drop=F]
      osi<-offset[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i])]
    }

    if(family == 'Bernoulli'){
      pij<-as.vector(boot::inv.logit(osi + xi%*%bta))
      u.ui<-pbinom(yi,size=1,prob=pij)
      l.ui<-pbinom(yi-1,size=1,prob=pij)
    }

    if(family == 'Poisson'){
      phi.theta<-as.vector(exp(osi + xi%*%bta))
      u.ui<-ppois(yi,phi.theta)
      l.ui<-ppois(yi-1,phi.theta)
    }

    if(family == 'Nbinom'){
      mui<-as.vector(exp(osi + xi%*%bta))
      sizei<-delta
      probi<-sizei/(sizei+mui)
      u.ui<-pnbinom(yi,size=sizei,prob=probi)
      l.ui<-pnbinom(yi-1,size=sizei,prob=probi)
    }

    zi.up<-qnorm(u.ui)
    zi.lo<-qnorm(l.ui)

    Ti<-phi<-matrix(0,m[i],m[i])
    if(m[i]>1){
      Ti[1,1]<-1
      for(j in 2:m[i]){
        for(k in 1:(j-1)){
          phi[j,k]<-pi/2-atan(W[i,j,k,]%*%gma)
        }
        Ti[j,j]<-prod(sin(phi[j,1:(j-1)]))
        Ti[j,1]<-cos(phi[j,1])
        if(j>2){
          for(l in 2:(j-1))
            Ti[j,l]<-cos(phi[j,l])*prod(sin(phi[j,1:(l-1)]))
        }
      }
      Ri<-Ti%*%t(Ti)
      temp<-0
      Rkl<-diag(2)
      for(l in 2:m[i])
        for(k in 1:(l-1)){
          Rkl[1,2]<-Rkl[2,1]<-Ri[k,l]
          prob<-mnormt::biv.nt.prob(df=Inf,lower=zi.lo[c(k,l)],upper=zi.up[c(k,l)],mean=c(0,0),S=Rkl)
          if(prob<=0) return(Inf)
          temp<-temp+log(prob)
        }
    }else{
      if(family == 'Bernoulli'){
        temp<-log(dbinom(yi,size=1,prob=pij))
      }

      if(family == 'Poisson'){
        temp<-log(dpois(yi,phi.theta))
      }

      if(family == 'Nbinom'){
        temp<-log(dnbinom(yi,size=sizei,prob=probi))
      }
    }

    lik<-lik+temp
  }
  return(-lik)
}



##--------------------------------------------------------------------------------------------------
##gmalik
gmalik <- function(gma, bta.delta, m, Y, X, W, family, offset)
{
  #Y: n x max(m) array
  #X: N x p covariates matrix
  #W: n x max(m) x (max(m)-1) x length(gamma) array
  #m: measurement vector of lengt n

  n<-length(m)
  p<-ncol(X)
  bta<-bta.delta[1:p]

  if(family == 'Nbinom'){
    bta<-bta.delta[1:p]
    delta<-bta.delta[p+1]
  }

  lik<-0
  for(i in 1:n){
    yi<-Y[i,1:m[i]]
    if(i==1){
      xi<-X[1:m[1],,drop=F]
      osi<-offset[1:m[1]]
    }else{
      xi<-X[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i]),,drop=F]
      osi<-offset[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i])]
    }

    if(family == 'Bernoulli'){
      pij<-as.vector(boot::inv.logit(osi + xi%*%bta))
      u.ui<-pbinom(yi,size=1,prob=pij)
      l.ui<-pbinom(yi-1,size=1,prob=pij)
    }

    if(family == 'Poisson'){
      phi.theta<-as.vector(exp(osi + xi%*%bta))
      u.ui<-ppois(yi,phi.theta)
      l.ui<-ppois(yi-1,phi.theta)
    }

    if(family == 'Nbinom'){
      mui<-as.vector(exp(osi + xi%*%bta))
      sizei<-delta
      probi<-sizei/(sizei+mui)
      u.ui<-pnbinom(yi,size=sizei,prob=probi)
      l.ui<-pnbinom(yi-1,size=sizei,prob=probi)
    }

    zi.up<-qnorm(u.ui)
    zi.lo<-qnorm(l.ui)

    Ti<-phi<-matrix(0,m[i],m[i])
    if(m[i]>1){
      Ti[1,1]<-1
      for(j in 2:m[i]){
        for(k in 1:(j-1)){
          phi[j,k]<-pi/2-atan(W[i,j,k,]%*%gma)
        }
        Ti[j,j]<-prod(sin(phi[j,1:(j-1)]))
        Ti[j,1]<-cos(phi[j,1])
        if(j>2){
          for(l in 2:(j-1))
            Ti[j,l]<-cos(phi[j,l])*prod(sin(phi[j,1:(l-1)]))
        }
      }
      Ri<-Ti%*%t(Ti)
      temp<-0
      Rkl<-diag(2)
      for(l in 2:m[i])
        for(k in 1:(l-1)){
          Rkl[1,2]<-Rkl[2,1]<-Ri[k,l]
          prob<-mnormt::biv.nt.prob(df=Inf,lower=zi.lo[c(k,l)],upper=zi.up[c(k,l)],mean=c(0,0),S=Rkl)
          if(prob<=0) return(Inf)
          temp<-temp+log(prob)
        }
    }else{
      if(family == 'Bernoulli'){
        temp<-log(dbinom(yi,size=1,prob=pij))
      }

      if(family == 'Poisson'){
        temp<-log(dpois(yi,phi.theta))
      }

      if(family == 'Nbinom'){
        temp<-log(dnbinom(yi,size=sizei,prob=probi))
      }
    }

    lik<-lik+temp
  }
  return(-lik)
}




##----------------------------------------------------------------------------
##asycov
##var（theta_hat）
asycov <- function(m, Y, X, W, family, theta , offset)
{
  #Y: vector of size sum(m)
  #X: matrix of size sum(m)x lbta
  #W: array of size nx max(m)x max(m) x lgma

  if(family=="Nbinom"){
    Gmat_Inv <- NbiGmatinv(theta, m, Y, X, W, offset)
  }

  if(family=="Bernoulli"){
    Gmat_Inv <- BnGmatinv(theta, m, Y, X, W, offset)
  }

  if(family=="Poisson"){
    Gmat_Inv <- PoiGmatinv(theta, m, Y, X, W, offset)
  }

  return(Gmat_Inv)
}



###----------------------------------------------------------------------------
cat.f <- function(...){ cat(..., fill = TRUE)}



##-----------------------------------------------------------------------------
##---------------------------------------------------------------------------
### LRT test: H0: theta0=theta0, theta=(theta0,theta1)
###---------------------------------------------------------------------------------------------
##id:the id of paremeter to test,in vector
##lrt.full: full model


#p.test
p.test <- function(id, theta0, opt, std, n, p, family=family){
  theta = opt$par[id]
  tval = (theta-theta0)/std[id]
  pval = 1 - pt(tval, df = n-p)
  return(list(tval=tval, pval=pval))
}
