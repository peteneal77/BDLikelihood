#
# Generic MCMC code for an epidemic in progress
#
#
#

library(MASS) # MASS library used for sampling from multivariate normal distribution for the random walk Metropolis algorithm

# Load in functions for calculating log-likelihood

source("aGSEI_functions.R")
source("EGSEI_functions.R")
source("GSEI_functions.R")

# Random walk Metropolis algorithm

MCMC_IncGSE=function(t,st,flike,N,beta,gamma,burn,bnits,nits,sigma,nub=10,lamb=500/6,nug=10,lamg=100)
  # t - inter-removal times 
  # st - observation up to time st (if removal times are after st they are ignored)
  # flike - choice of likelihood functions
  # beta,gamma - initial parameter values
  # burn, bnits - burn is the number of burn-in periods each of length bnits 
  # nits - number of iterations for final MCMC run
  # sigma - sets initial (co)variance matrix for RWM with Sig = sigma*I
  # Gamma priors Gamma (nub,lamb) for beta and Gamma(nug,lamg) for gamma
{
  logl=flike(beta,gamma,t,st,N)+dgamma(beta,nub,lamb,log=TRUE)+dgamma(gamma,nug,lamg,log=TRUE) 
  # Initial calculation of log-posterior
  Sig=diag(rep(sigma^2,2)) # Set initial covariance matrix
  for(ii in 1:burn)
  {
    OUTX=matrix(0,ncol=2,nrow=bnits)
    for(ij in 1:bnits)
    {
      U=mvrnorm(1,rep(0,2),Sig)
      betanew=beta+U[1]
      gammanew=gamma+U[2]
      if(min(betanew,gammanew)>0) loglnew=flike(betanew,gammanew,t,st,N)+dgamma(betanew,nub,lamb,log=TRUE)+dgamma(gammanew,nug,lamg,log=TRUE)
      if(min(betanew,gammanew)<0) loglnew=-100000
      # If betanew or gammanew are negative the proposed move is rejected. 
      u=log(runif(1))
      if(u<(loglnew-logl))
      {
        beta=betanew
        gamma=gammanew
        logl=loglnew
      }
      OUTX[ij,]=c(beta,gamma)
      #      print(OUTX[ij,])
    }
    Svar=var(OUTX)
    Sig=(2.4^2/2)*(0.95*Svar+0.05*diag(diag(Svar))) # Updating of the proposal (co)variance matrix
    print(ii)
  }
  OUT=matrix(0,ncol=2,nrow=nits)
  for(ij in 1:nits)
  {
    U=mvrnorm(1,rep(0,2),Sig)
    betanew=beta+U[1]
    gammanew=gamma+U[2]
    if(min(betanew,gammanew)>0) loglnew=flike(betanew,gammanew,t,st,N)+dgamma(betanew,nub,lamb,log=TRUE)+dgamma(gammanew,nug,lamg,log=TRUE)
    if(min(betanew,gammanew)<0) loglnew=-100000
    u=log(runif(1))
    if(u<(loglnew-logl))
    {
      beta=betanew
      gamma=gammanew
      logl=loglnew
    }
    OUT[ij,]=c(beta,gamma)
  }
  OUT
}
