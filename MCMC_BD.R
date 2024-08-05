#
# BD MCMC algorithm
#
#
#

library(MASS) # Allows for simulation from multivariate normal distribution

source("BD_Exact_Like.R") # Exact log-likelihood
source("BD_Approx_Like.R") # Approximate log-likelihood
source("BD_Hybrid_Like_Fix.R") # Hybrid log-likelihood - with fixed change point from exact to approximate likelihood

# 

source("BD_para.R") # Converts parameters theta to alpha, mu, d

source("BD_prior.R") # Computes prior density for parameters theta

# General likelihood function call which has input alpha, mu, d
# Parameter function to change parameters into alpha, mu, d
#


MCMC_BD=function(t,flike,fpara,theta,Nrun,sigtuna,theta_type,theta_prior,theta_fix)
  # t - inter-removal times 
  # flike - choice of likelihood functions
  # fpara - choice of parameter function theta -> alpha, mu, d
  # theta - initial parameter values
  # Nrun - Number of iterations per run (all but last are burn-in)   
  # nits - number of iterations for final MCMC run
  # sigtuna - sets initial (co)variance matrix for RWM with Sig = sigtuna*I
  # theta_type defines type of prior on parameters theta 
  # (0 - real, 1 - positive, 2 - probability, 3 - integer on [0,N])
  # theta_prior defines the values of the prior distribution
  # - matrix (length(theta) * 3) - first entry type, second and third prior
  # - Types 1 - uniform (includes discrete uniform), 2 - gamma, 3 - beta, 4 - normal
  # theta_fix (0 - fixed, 1 - updated)
{
  thetanew=theta # Set new values equal to theta - fixed values will not update
  
  nl=length(Nrun) # Number of loops of MCMC run
  # After each loop the proposal variance for the RWM algorithm is updated.
  
  np=sum(theta_fix) # Number of parameters which are updated
  npa=length(theta) # Total number of parameters

  NN=length(t)+1  # Total number of detected deaths
  
  Pa=fpara(theta,NN) # Generates alpha, mu, d from the parameters theta
  alpha=Pa$alpha
  mu=Pa$mu
  d=Pa$d
  logl=flike(alpha,mu,d,t) #  Initial calculation of log-likelihood 
  logp=priorden(theta,theta_prior)# Initial calculation of log-prior
  ss=rep(sigtuna^2,np) # Set initial variances for RWM
  for(i in 1:npa)
  {
    if((theta_fix[i]==1)&(theta_type[i]==3)) ss[i]=1 # For discrete parameters set variance equal to 1
  }
  Sig=diag(ss) # Set initial covariance matrix
  for(ii in 1:nl)
  {
    acc=0 # Used to calculate acceptance rate
    OUTX=matrix(0,ncol=npa,nrow=Nrun[ii]) # Stores each MCMC iteration
    for(ij in 1:Nrun[ii])
    {
      U=mvrnorm(1,rep(0,np),Sig) # multivariate normal proposal for updated parameters
      thetanew[theta_fix==1]=theta[theta_fix==1]+U # Update parameters

      thetanew[theta_type==3]=round(thetanew[theta_type==3],0) # Round discrete parameters to integer values
      
#      print(thetanew)
      
      logpnew=priorden(thetanew,theta_prior) # Compute log-prior
      
      if(logpnew>-Inf) # If parameters consistent with prior
      {
        Pa=fpara(thetanew,N)  # Compute updated alpha, mu, d values
        alphanew=Pa$alpha
        munew=Pa$mu
        dnew=Pa$d
      
        loglnew=flike(alphanew,munew,dnew,t) # Compute log-likelihood for thetanew
     
         u=log(runif(1))
         if(u<(loglnew+logpnew-logl-logp)) # Acceptance decision
         {
           acc=acc+1 # Counts number of accepted values
           theta=thetanew # update theta
           logl=loglnew   # update log-likelihood
           logp=logpnew   # update log-prior
         }
      }   
      OUTX[ij,]=theta # Store theta values
    }
    Svar=var(OUTX[theta_fix==1,theta_fix==1])
    Sig=(2.4^2/np)*(0.95*Svar+0.05*diag(diag(Svar))) # Updating of the proposal (co)variance matrix
    print(ii) # Allows us to keep track of MCMC run
    print(acc/Nrun[ii]) # Proportion of accepted values in run
  }
  OUTX
}
