#
# Approximate likelihood for time-inhomogeneous Birth-Death process using Gaussian approximations
# - Piecewise constant parameters between death times
# - Partial detection of deaths (allowed for)

# alpha_k, mu_k, d_k - birth and death rates and detection probability between deaths (k-1) and k
#

# xij denotes the j^th derivative of \xi given in (3.26) and is used to find the mean and variance of the
# Gaussian approximation in (3.27) and (3.28)

xi0=function(kappa,sigma,theta)
{
  exp(theta*kappa+0.5*theta^2*sigma)
}

xi1=function(kappa,sigma,theta)
{
  A=exp(theta*kappa+0.5*theta^2*sigma)
  (kappa+theta*sigma)*A
}

xi2=function(kappa,sigma,theta)
{
  A=exp(theta*kappa+0.5*theta^2*sigma)
  (sigma+(kappa+theta*sigma)^2)*A
}

xi3=function(kappa,sigma,theta)
{
  A=exp(theta*kappa+0.5*theta^2*sigma)
  (3*sigma+(kappa+theta*sigma)^2)*(kappa+theta*sigma)*A
}



#
# alpha, mu 
# Approximate likelihood calculation - Birth-Death process is ongoing
#

alike_BD=function(alpha,mu,d,t,taut)
# Parameters: alpha - birth, mu - death, d - detection prob
# t - inter-arrival times of detected deaths
# taut - time since the last death
# K - deaths means t is of length K-1 and alpha, mu, d are of length K+1  
{
  tt=c(0,t,taut) # Inter-arrival times + time since last death + 0 at the front as place holder
  K=length(t)+1 #  Detected number of deaths
  p=alpha/(alpha+mu) # Birth probability p_k - between detected deaths k-1 and k 
  q=1-p # Death probability p_k - between detected deaths k-1 and k
  u=sqrt(1-4*p*q*(1-d))  # u_k
  lambda=(1+u-2*p)/(1+u)  # lambda_k
  nu=(1-u)/(2*p)  # nu_k
  phi=exp(-(alpha+mu)*u*tt)  # phi_k - do not use phi_1
  psi=(1-lambda)*(1-phi)/(1-nu*(1-lambda)*phi)  # psi_k - do not use psi_1
  pi=lambda[1] # Initial value of pi
  for(k in 2:(K+1))
  {
    pi[k]=(lambda[k]*(1-nu[k]*(1-pi[k-1]))-(1-nu[k])*(lambda[k]-pi[k-1])*phi[k])/
      (1-nu[k]*(1-pi[k-1])+nu[k]*(lambda[k]-pi[k-1])*phi[k])
  }
  # Sets value of pi
  r=1 # Initial r not used
  r[2:(K+1)]=pi[1:K]*(lambda[2:(K+1)]+(1-lambda[2:(K+1)])*(1-nu[2:(K+1)])*phi[2:(K+1)])/
    (lambda[2:(K+1)]*(1-nu[2:(K+1)]*(1-pi[1:K]))-(1-nu[2:(K+1)])*(lambda[2:(K+1)]-pi[1:K])*phi[2:(K+1)])
  # Sets value of r
  g=0 # Not used but g[2:(K+1)] is used
  g[2:(K+1)]=pi[1:K]*lambda[2:(K+1)]*((1-nu[2:(K+1)])*(lambda[2:(K+1)]-pi[1:K])
              +(1-lambda[2:(K+1)])*(1-nu[2:(K+1)])*(1-nu[2:(K+1)]*(1-pi[1:K])))*
              (alpha[2:(K+1)]+mu[2:(K+1)])*u[2:(K+1)]*phi[2:(K+1)]/
              (lambda[2:(K+1)]*(1-nu[2:(K+1)]*(1-pi[1:K]))-(1-nu[2:(K+1)])*(lambda[2:(K+1)]-pi[1:K])*phi[2:(K+1)])^2
  # Calculates g
  h=(1-pi-psi)/((1-pi)*(1-psi)) # Computes h
  logl=log(g[2]) # initial value of log-likelihood
  
  eta=2 # initial values for eta and sigma corresponding to R_2 = 2
  sigma=0
  
  for(k in 3:K)
  {
    theta=log(r[k])  # Calculate log r_k
    c1=xi1(eta,sigma,theta)
    c2=xi2(eta,sigma,theta)
    c3=xi3(eta,sigma,theta)
    # Compute \xi_k^(i)  i=1,2,3
    
    logl=logl+log(g[k])+log(c1)-log(r[k])  # Updates log-likelihood with contribution for removal k
    
    eta=(2*c1+h[k]*(c2-c1))/c1 # Update mean
    sigma=((4*c1+5*h[k]*(c2-c1)+h[k]^2*(c3-3*c2+2*c1))/c1)-eta^2  # Update variance
  }
  c0=xi0(eta,sigma,log(r[K+1])) 
  logl=logl+log(c0) # Log-likelihood for nothing happening since last removal.
  logl
}



