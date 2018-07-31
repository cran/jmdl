##calculate Gmat.inv for three model(Nbinom, Bernoulli, Poisson)

##for Nbinom
NbiGmatinv <- function(theta, m, Y, X, W, offset)
{
  n<-length(m)
  lth<-length(theta)
  Kmat<-matrix(0,lth,lth)
  Jmat<-matrix(0,lth,lth)

  p<-ncol(X)
  bta<-theta[1:p]
  delta<-theta[p+1]
  gma<-theta[-c(1:(p+1))]

  for(i in 1:n){
    if(i==1){
      xi<-X[1:m[1],,drop=F]
      osi <- offset[1:m[1]]
    }else{
      xi<-X[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i]),,drop=F]
      osi<-offset[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i])]
    }
    yi<-Y[i,1:m[i]]
    Ti<-phi<-diag(1,m[i],m[i])
    if(m[i]>1){
      for(j in 2:m[i]){
        for(k in 1:(j-1))
          phi[j,k]<-pi/2-atan(W[i,j,k,]%*%gma)

        Ti[j,j]<-prod(sin(phi[j,1:(j-1)]))
        Ti[j,1]<-cos(phi[j,1])
        if(j>2){
          for(l in 2:(j-1))
            Ti[j,l]<-cos(phi[j,l])*prod(sin(phi[j,1:(l-1)]))
        }
      }
    }
    Ri<-Ti%*%t(Ti)
    Sni<-rep(0,lth)
    if(m[i]==1){
      llik.dev<-c(nb.deriv.pmf(bta,delta,yi,xi,osi,'Nbinom'),rep(0,length(gma)))
      Sni<-Sni+llik.dev
      Kmat<-Kmat+llik.dev%*%t(llik.dev)
    }else{
      wi<-W[i,,,]
      for(k in 2:m[i])
        for(j in 1:(k-1)){
          llik.dev<-nb.deriv.lijk(j,k,bta,delta,gma,Ri[j,k],yi,xi,osi,wi,phi,Ti,'Nbinom')
          Sni<-Sni+llik.dev
          Kmat<-Kmat+llik.dev%*%t(llik.dev)
        }
    }
    Jmat<-Jmat+Sni%*%t(Sni)
  }
  Kmat.inv<-MASS::ginv(Kmat)
  Gmat.inv<-(Kmat.inv%*%Jmat)%*%Kmat.inv
  return(list(Gmat.inv=Gmat.inv,Kmat.inv=Kmat.inv,Jmat=Jmat))
}



##Gmat.inv
##for Bernoulli
BnGmatinv <- function(theta, m, Y, X, W, offset)
{
  n<-length(m)
  lth<-length(theta)
  Kmat<-matrix(0,lth,lth)
  Jmat<-matrix(0,lth,lth)

  p<-ncol(X)
  bta<-theta[1:p]
  gma<-theta[-c(1:p)]

  for(i in 1:n){
    if(i==1){
      xi<-X[1:m[1],,drop=F]
      osi <- offset[1:m[1]]
    }else{
      xi<-X[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i]),,drop=F]
      osi<-offset[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i])]
    }
    yi<-Y[i,1:m[i]]
    Ti<-phi<-matrix(0,m[i],m[i])
    Ti[1,1]<-1
    if(m[i]>1){
      for(j in 2:m[i]){
        for(k in 1:(j-1))
          phi[j,k]<-pi/2-atan(W[i,j,k,]%*%gma)
        Ti[j,j]<-prod(sin(phi[j,1:(j-1)]))
        Ti[j,1]<-cos(phi[j,1])
        if(j>2){
          for(l in 2:(j-1))
            Ti[j,l]<-cos(phi[j,l])*prod(sin(phi[j,1:(l-1)]))
        }
      }
    }
    Ri<-Ti%*%t(Ti)
    Sni<-rep(0,lth)
    if(m[i]==1){
      llik.dev<-c(deriv.pmf(bta,yi,xi,osi,"Bernoulli"),rep(0,length(gma)))
      Sni<-Sni+llik.dev
      Kmat<-Kmat+llik.dev%*%t(llik.dev)
    }else{
      wi<-W[i,,,]
      for(k in 2:m[i])
        for(j in 1:(k-1)){
          llik.dev<-deriv.lijk(j,k,bta,gma,Ri[j,k],yi,xi,osi,wi,phi,Ti,"Bernoulli")
          if(any(is.na(llik.dev))) cat(c(i,j,k),"\n")
          Sni<-Sni+llik.dev
          Kmat<-Kmat+llik.dev%*%t(llik.dev)
        }
    }
    Jmat<-Jmat+Sni%*%t(Sni)
  }
  Kmat.inv<-MASS::ginv(Kmat)
  Gmat.inv<-(Kmat.inv%*%Jmat)%*%Kmat.inv
  return(list(Gmat.inv=Gmat.inv,Kmat.inv=Kmat.inv,Jmat=Jmat))
}



##Gmat.inv
##for Poisson
PoiGmatinv <- function(theta, m, Y, X, W, offset)
{
  n<-length(m)
  lth<-length(theta)
  Kmat<-matrix(0,lth,lth)
  Jmat<-matrix(0,lth,lth)
  p<-ncol(X)
  bta<-theta[1:p]
  gma<-theta[-c(1:p)]

  for(i in 1:n){
    if(i==1){
      xi<-X[1:m[1],,drop=F]
      osi<-offset[1:m[1]]
    }else{
      xi<-X[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i]),,drop=F]
      osi<-offset[(sum(m[1:(i-1)])+1):(sum(m[1:(i-1)])+m[i])]
    }
    yi<-Y[i,1:m[i]]
    Ti<-phi<-diag(1,m[i],m[i])
    if(m[i]>1){
      for(j in 2:m[i]){
        for(k in 1:(j-1))
          phi[j,k]<-pi/2-atan(W[i,j,k,]%*%gma)

        Ti[j,j]<-prod(sin(phi[j,1:(j-1)]))
        Ti[j,1]<-cos(phi[j,1])
        if(j>2){
          for(l in 2:(j-1))
            Ti[j,l]<-cos(phi[j,l])*prod(sin(phi[j,1:(l-1)]))
        }
      }
    }
    Ri<-Ti%*%t(Ti)
    Sni<-rep(0,lth)
    if(m[i]==1){
      llik.dev<-c(deriv.pmf(bta,yi,xi,osi,'Poisson'),rep(0,length(gma)))
      Sni<-Sni+llik.dev
      Kmat<-Kmat+llik.dev%*%t(llik.dev)
    }else{
      wi<-W[i,,,]
      for(k in 2:m[i])
        for(j in 1:(k-1)){
          llik.dev<-deriv.lijk(j,k,bta,gma,Ri[j,k],yi,xi,osi,wi,phi,Ti,"Poisson")
          Sni<-Sni+llik.dev
          Kmat<-Kmat+llik.dev%*%t(llik.dev)
        }
    }
    Jmat<-Jmat+Sni%*%t(Sni)
  }
  Kmat.inv<-MASS::ginv(Kmat)
  Gmat.inv<-(Kmat.inv%*%Jmat)%*%Kmat.inv
  return(list(Gmat.inv=Gmat.inv,Kmat.inv=Kmat.inv,Jmat=Jmat))
}
