#
# GSE - approximation using Gaussian
# - Constructed and run Abikiliki data
#  Set d=1
#

# Time-inhomogeneous birth-death process 
# - Piecewise constant parameters between death times
# - Partial detection of deaths (allowed for)

# alpha_k, mu_k, - birth and death rates and detection probability
#
# Returns log-likelihood

# Calculations of xi_k (k=0,1,2,3)
#
#

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
# beta, gamma 
#
#

alike_GSEI=function(beta,gamma,t,st,N)
  # Vectors of parameters beta, gamma
  # t - inter-arrival times between events
  # st - stopping time
  # N - population size
{
  s=t
  for(i in 1:length(t)) s[i]=sum(t[1:i]) # time of the i^{th} removal
  t=t[s<=st]
  s=s[s<=st]
  tau=st-max(s) # time since the last removal
  
  # Above sets up inference based on removal times upto time st.  
  
  alpha=beta*(N-1)/N # Initial birth rate equal to infection rate * proportion of susceptibles
  k=length(t)+1 # Total number of deaths 
  p=alpha/(alpha+gamma)  # Probability an event is a birth
  q=1-p        # Probability an event is a death
  pi=q       # Set pi = lambda
  
  # First infection
  
  infect=(1-pi)/pi # Estimated number infected
  alpha=max(beta*(N-(i-1)-infect)/N,0) # Estimated proportion susceptible
  p=alpha/(alpha+gamma)  # Probability an event is a birth
  q=1-p        # Probability an event is a death
  phi=exp(-(alpha+gamma)*t[1])
  psi=(1-q)*(1-phi)
  g=pi*q*(1-pi)*(alpha+gamma)*phi/(q-(q-pi)*phi)^2
  r=pi*(q + (1-q )*phi)/(q -(q -pi)*phi)
  # pi update
  pi=q-(q-pi)*phi
  h=(1-pi-psi)/((1-pi)*(1-psi)) # Uses the updated pi
  logl=log(g)  # 
           
  kappa=2 # initial values for kappa and sigma corresponding to R_2 = 2
  sigma=0
           
  for(i in 3:k)
  {
    infect=kappa*(1-pi)/pi # Estimated number infected - kappa = E[R_k] (approx)
#    alpha=beta*(N-(i-1)-infect)/N # i-1 individuals have been removed
    alpha=max(beta*(N-(i-1)-infect)/N,0) # alpha capped to be non-negative
    p=alpha/(alpha+gamma)  # Probability an event is a birth
    q=1-p        # Probability an event is a death
    phi=exp(-(alpha+gamma)*t[i-1])
    psi=(1-q )*(1-phi)
    g=pi*q*(1-pi)*(alpha+gamma)*phi/(q-(q-pi)*phi)^2
    r=pi*(q + (1-q )*phi)/(q -(q -pi)*phi)
                      # pi update
    pi=q-(q-pi)*phi
    h=(1-pi-psi)/((1-pi)*(1-psi)) # Uses the updated pi
                      
    theta=log(r)
    c0=xi0(kappa,sigma,theta)
    c1=xi1(kappa,sigma,theta)
    c2=xi2(kappa,sigma,theta)
    c3=xi3(kappa,sigma,theta)
                      
    logl=logl+log(g)+log(c1)-log(r)
                      
    kappa=(2*c1+h*(c2-c1))/c1
    sigma=((4*c1+5*h*(c2-c1)+h^2*(c3-3*c2+2*c1))/c1)-kappa^2
  }  
  # After last observed removal
  
  infect=kappa*(1-pi)/pi # Estimated number infected
  alpha=max(beta*(N-(i-1)-infect)/N,0)
  p=alpha/(alpha+gamma)
  q=1-p
  phi=exp(-(alpha+gamma)*tau)
  r=pi*(q+(1-q)*phi)/(q-(q-pi)*phi)
  theta=log(r)
  logl=logl+theta*kappa+0.5*theta^2*sigma
  logl
}



