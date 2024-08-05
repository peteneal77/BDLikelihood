#
# Birth-death process likelihood
# 

# Hybrid likelihood switch from exact to approximate after the L^th detected death.
# 
# 01.03.24 - Code adapted to calculate the likelihood for each L (avoids recalculating full likelihood each time).

# alpha - birth rates
# mu - death rates
# d - probability of detecting a death

# Compute likelihood given inter-arrival times using t.

like_BD_hybrid_ext=function(alpha,mu,d,t)
  # Vectors of birth, death rates and probability
  # t - inter-arrival times between events
{
  
  s=t
  N=length(t)+1 # Number of removals as t =(t_2, ..., t_N)
  
  Lv=rep(0,N) # Stores calculation of log-likelihood
  
  for(i in 1:N) s[i]=sum(t[1:i]) # time of the i^{th} removal
  
  p=alpha/(mu+alpha) # probability event is a birth
  q=1-p # probability event is a death
  
  # Quantities from (8) in paper.
  
  u=sqrt(1-4*p*q*(1-d))
  lambda=(1+u-2*p)/(1+u)
  nu=(1-u)/(2*p)
  phik=exp(-(alpha+mu)*u*c(0,t)) # Added 0 to create phi_k
  psi=(1-lambda)*(1-phik)/(1-nu*(1-lambda)*phik)
  
  # logL = log Likelihood  
  
  logL=0  
  
  # pi, r and g are needed in the main loop   
  pi=lambda[1]
  
  # Log-likelihood k=2
  
  k=2
  gnum=pi*lambda[k]*(1-nu[k])*(lambda[k]-pi+(1-lambda[k])*(1-nu[k]*(1-pi)))*(alpha[k]+mu[k])*u[k]*phik[k]
  gden=lambda[k]*(1-nu[k]*(1-pi))-(1-nu[k])*(lambda[k]-pi)*phik[k]
  
  logL=log(gnum)-2*log(gden)
  
  # Update pi
  pi=(lambda[k]*(1-nu[k]*(1-pi))-(1-nu[k])*(lambda[k]-pi)*phik[k])/(1-nu[k]*(1-pi)+nu[k]*(lambda[k]-pi)*phik[k])
  
  B=1 # Distribution of R_2, a point mass at R_2 =2
  
  eta=2 # Distribution of R_2, a point mass at R_2 =2
  sigma=0
  
  logE=logL
  
  pi1=pi
  for(j in 3:N)
  {
    r1=pi1*(lambda[j]+(1-lambda[j])*(1-nu[j])*phik[j])/(lambda[j]*(1-nu[j]*(1-pi1))-(1-nu[j])*(lambda[j]-pi1)*phik[j])
    
    gnum1=pi1*lambda[j]*(1-nu[j])*(lambda[j]-pi1+(1-lambda[j])*(1-nu[j]*(1-pi1)))*(alpha[j]+mu[j])*u[j]*phik[j]
    gden1=lambda[j]*(1-nu[j]*(1-pi1))-(1-nu[j])*(lambda[j]-pi1)*phik[j]
    
    theta=log(r1)
    
    c1=xi1(eta,sigma,theta)
    
    
    logE=logE+log(c1)-log(r1)+log(gnum1)-2*log(gden1)
    
    # Update pi, h eta and sigma   
    pi1=(lambda[j]*(1-nu[j]*(1-pi1))-(1-nu[j])*(lambda[j]-pi1)*phik[j])/(1-nu[j]*(1-pi1)+nu[j]*(lambda[j]-pi1)*phik[j])
    
    h1=(1-pi1-psi[j])/((1-pi1)*(1-psi[j]))
    
    eta0=eta # current
    sigma0=sigma 
    
    eta=2+h1*(eta0+theta*sigma0+sigma0/(eta0+theta*sigma0)-1)
    sigma=h1*(1-h1)*(eta0+theta*sigma0+sigma0/(eta0+theta*sigma0)-1)+h1^2*sigma0*(1-sigma0/(eta0+theta*sigma0)^2)
  }
  
  Lv[2]=logE
  
  for(k in 3:N)
  {
    r=pi*(lambda[k]+(1-lambda[k])*(1-nu[k])*phik[k])/(lambda[k]*(1-nu[k]*(1-pi))-(1-nu[k])*(lambda[k]-pi)*phik[k])
    ss=seq(2,k-1) # Possible values for R_{k-1}
    Er=sum(B*ss*r^(ss-1)) 
    
    gnum=pi*lambda[k]*(1-nu[k])*(lambda[k]-pi+(1-lambda[k])*(1-nu[k]*(1-pi)))*(alpha[k]+mu[k])*u[k]*phik[k]
    gden=lambda[k]*(1-nu[k]*(1-pi))-(1-nu[k])*(lambda[k]-pi)*phik[k]
    
    logL=logL+log(Er)+log(gnum)-2*log(gden) # Exact log-likelihood
    
    # Update pi and B    
    pi=(lambda[k]*(1-nu[k]*(1-pi))-(1-nu[k])*(lambda[k]-pi)*phik[k])/(1-nu[k]*(1-pi)+nu[k]*(lambda[k]-pi)*phik[k])
    h=(1-pi-psi[k])/((1-pi)*(1-psi[k]))
    M=Mmat(h,r,k-1)
    B=B%*%M
    B=B/sum(B)
    
    BJ=seq(2,k,1)
    eta=sum(B*BJ)
    sigma=sum(B*BJ^2)-eta^2
    
    logE=logL
    pi1=pi
    if(k<N)
    {  
    for(j in (k+1):N)
    {
      r1=pi1*(lambda[j]+(1-lambda[j])*(1-nu[j])*phik[j])/(lambda[j]*(1-nu[j]*(1-pi1))-(1-nu[j])*(lambda[j]-pi1)*phik[j])
      
      gnum1=pi1*lambda[j]*(1-nu[j])*(lambda[j]-pi1+(1-lambda[j])*(1-nu[j]*(1-pi1)))*(alpha[j]+mu[j])*u[j]*phik[j]
      gden1=lambda[j]*(1-nu[j]*(1-pi1))-(1-nu[j])*(lambda[j]-pi1)*phik[j]
      
      theta=log(r1)
      
      c1=xi1(eta,sigma,theta)
      
      
      logE=logE+log(c1)-log(r1)+log(gnum1)-2*log(gden1)
      
      # Update pi, h eta and sigma   
      pi1=(lambda[j]*(1-nu[j]*(1-pi1))-(1-nu[j])*(lambda[j]-pi1)*phik[j])/(1-nu[j]*(1-pi1)+nu[j]*(lambda[j]-pi1)*phik[j])
      
      h1=(1-pi1-psi[j])/((1-pi1)*(1-psi[j]))
      
      eta0=eta # current
      sigma0=sigma 
      
      eta=2+h1*(eta0+theta*sigma0+sigma0/(eta0+theta*sigma0)-1)
      sigma=h1*(1-h1)*(eta0+theta*sigma0+sigma0/(eta0+theta*sigma0)-1)+h1^2*sigma0*(1-sigma0/(eta0+theta*sigma0)^2)
    }
    }
    Lv[k]=logE
  }
  

  #  print(eta)
  #  print(sigma)
  #  print(logL)
  Lv
}  

