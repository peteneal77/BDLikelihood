#
# Generic MCMC code for an epidemic in progress
#
# Stores the predicted number infected at the end of the observation period (last removal before end)
# Code modified from inter-removal times to daily count data.
# All deaths are detected

library(MASS) # MASS library used for sampling from multivariate normal distribution for the random walk Metropolis algorithm

# Load in functions for calculating log-likelihood

source("aGSEI_functions.R")
source("EGSEI_functions.R")
source("GSEI_functions.R")
source("daily_convert.R")

# Random walk Metropolis algorithm

MCMC_daily_IncGSE=function(rem,st,flike,N,beta,gamma,burn,bnits,nits,sigma,Nup,fixp,
                      nub=10,lamb=500/6,nug=10,lamg=100)
  # rem - Number of removal times per day
  # st - observation up to time st (if removal times are after st they are ignored)
  # flike - choice of likelihood functions
  # beta,gamma - initial parameter values
  # burn, bnits - burn is the number of burn-in periods each of length bnits 
  # nits - number of iterations for final MCMC run
  # sigma - sets initial (co)variance matrix for RWM with Sig = sigma*I
  # Gamma priors Gamma (nub,lamb) for beta and Gamma(nug,lamg) for gamma
{
  stx=length(rem)
  TDcase=rep(seq(1,stx),rem)
  
  # Initialise removal times
  DA=DailyCon(rem,0)
  t=DA$ainter
  tim=DA$acu
  
  logl=flike(beta,gamma,t,st,N)+dgamma(beta,nub,lamb,log=TRUE)+dgamma(gamma,nug,lamg,log=TRUE)
  # Initial calculation of log-posterior
  Sig=diag(rep(sigma^2,2)) # Set initial covariance matrix
  for(ii in 1:burn)
  {
    OUTX=matrix(0,ncol=2,nrow=bnits)
    upx=0
    for(ij in 1:bnits)
    {
      # Update removal times
      DAnew=InterUp(tim,TDcase,Nup)
      tnew=DAnew$ainter
      timnew=DAnew$acu
      loglnew=flike(beta,gamma,tnew,st,N)+dgamma(beta,nub,lamb,log=TRUE)+dgamma(gamma,nug,lamg,log=TRUE)

      u=log(runif(1))
      if(u<(loglnew-logl))
      {
        t=tnew
        tim=timnew
        logl=loglnew
        upx=upx+1
      }
      
      # Update parameters
      U=mvrnorm(1,rep(0,2),Sig)
      U=U*fixp
      betanew=beta+U[1]
      gammanew=gamma+U[2]
      if(min(c(betanew,gammanew))>0) 
      {
        loglnew=flike(betanew,gammanew,t,st,N)+dgamma(betanew,nub,lamb,log=TRUE)+dgamma(gammanew,nug,lamg,log=TRUE) 
      }
      if(min(c(betanew,gammanew))<0) loglnew=-100000
      u=log(runif(1))
      if(u<(loglnew-logl))
      {
        beta=betanew
        gamma=gammanew
        logl=loglnew
      }
      OUTX[ij,]=c(beta,gamma)
    }
    Svar=var(OUTX[,1:2])
    Sig=(2.4^2/2)*(0.95*Svar+0.05*diag(diag(Svar))) # Updating of the proposal (co)variance matrix
    Sig=Sig+diag((1-fixp)*0.01) # This allows for none problematic cases if a parameter is fixed
    print(ii)
    print(upx)
  }
  OUT=matrix(0,ncol=2,nrow=nits)
  upx=0
  for(ij in 1:nits)
  {
    # Update removal times
    DAnew=InterUp(tim,TDcase,Nup)
    tnew=DAnew$ainter
    timnew=DAnew$acu
    loglnew=flike(beta,gamma,tnew,st,N)+dgamma(beta,nub,lamb,log=TRUE)+dgamma(gamma,nug,lamg,log=TRUE) 
    u=log(runif(1))
    if(u<(loglnew-logl))
    {
      t=tnew
      tim=timnew
      logl=loglnew
      upx=upx+1
    }
    
    
    # Update parameters
    U=mvrnorm(1,rep(0,2),Sig)
    U=U*fixp
    betanew=beta+U[1]
    gammanew=gamma+U[2]
    if(min(c(betanew,gammanew))>0)  
    {
      loglnew=flike(betanew,gammanew,t,st,N)+dgamma(betanew,nub,lamb,log=TRUE)+dgamma(gammanew,nug,lamg,log=TRUE)
    }
    if(min(c(betanew,gammanew))<0) loglnew=-100000
    u=log(runif(1))
    if(u<(loglnew-logl))
    {
      beta=betanew
      gamma=gammanew
      logl=loglnew
    }
    OUT[ij,]=c(beta,gamma)
  }
  print(upx)
  OUT
}
