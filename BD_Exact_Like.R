#
# Birth-death process likelihood
# 

# Function to calculate the matrix M
# Note that this r * matrix given in paper.

Mmat=function(h,r,k)
{
  S=seq(0,k-1)
  M=matrix(0,nrow=k-1,ncol=k)
  for(i in 1:(k-1)) M[i,]=(i+1)*dbinom(S,i,h)*r^i
  M
}

# alpha - birth rates
# mu - death rates
# d - probability of detecting a death

# Compute likelihood given inter-arrival times using t.

like_BD_exact=function(alpha,mu,d,t)
# Vectors of birth, death rates and probability
# t - inter-arrival times between events
{
  s=t
  N=length(t)+1 # Number of removals as t =(t_2, ..., t_N)
  
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
  
  for(k in 3:N)
  {
    r=pi*(lambda[k]+(1-lambda[k])*(1-nu[k])*phik[k])/(lambda[k]*(1-nu[k]*(1-pi))-(1-nu[k])*(lambda[k]-pi)*phik[k])
    ss=seq(2,k-1) # Possible values for R_{k-1}
    Er=sum(B*ss*r^(ss-1)) 
    
    gnum=pi*lambda[k]*(1-nu[k])*(lambda[k]-pi+(1-lambda[k])*(1-nu[k]*(1-pi)))*(alpha[k]+mu[k])*u[k]*phik[k]
    gden=lambda[k]*(1-nu[k]*(1-pi))-(1-nu[k])*(lambda[k]-pi)*phik[k]
    
    logL=logL+log(Er)+log(gnum)-2*log(gden)
    
# Update pi and B    
    pi=(lambda[k]*(1-nu[k]*(1-pi))-(1-nu[k])*(lambda[k]-pi)*phik[k])/(1-nu[k]*(1-pi)+nu[k]*(lambda[k]-pi)*phik[k])
    h=(1-pi-psi[k])/((1-pi)*(1-psi[k]))
    M=Mmat(h,r,k-1)
    B=B%*%M
    B=B/sum(B)
  }
#  print(sum(B*seq(2,N)))   - Expected value of R_N
#  print(sum(B*seq(2,N)^2)-sum(B*seq(2,N))^2) - Variance of R_N
#  print(logL) - logL log likelihood
#  print(B)   - Distribution of R_N
  logL
}  
  

 