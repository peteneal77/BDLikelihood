#
# GSE - birth-death process approximation using ideas from Ball and Neal (2023)
# - Constructed and run Abikiliki data
# - Code simplified for case d=1 - all removals detected
#

# Time-inhomogeneous birth-death process 
# - Piece-wise constant parameters between death times
# - Partial detection of deaths (allowed for)

# alpha_k, mu_k - birth and death rates 
#
#

# Function to calculate the matrix M

Mmat=function(h,r,k)
{
  S=seq(0,k-1)
  M=matrix(0,nrow=k-1,ncol=k)
  for(i in 1:(k-1)) M[i,]=(i+1)*dbinom(S,i,h)*r^i
  M
}

#
# beta, gamma 
#
#

like_GSEI=function(beta,gamma,t,st,N)
  # Vectors of parameters beta and gamma 
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
  alpha=max(beta*(N-(i-1)-infect)/N,0) # sets proportion susceptible to have minimum = 0
  p=alpha/(alpha+gamma)  # Probability an event is a birth
  q=1-p        # Probability an event is a death
  phi=exp(-(alpha+gamma)*t[1])
  psi=(1-q)*(1-phi)
  g=pi*q*(1-pi)*(alpha+gamma)*phi/(q-(q-pi)*phi)^2
  r=pi*(q+ (1-q)*phi)/(q-(q-pi)*phi)
  # pi update
  pi=(q-(q-pi)*phi)
  h=(1-pi-psi)/((1-pi)*(1-psi)) # Uses the updated pi
  #  print(h)
  logl=log(g)  # 
  
  B=matrix(1,ncol=1,nrow=1)  # B matrix
  
  for(i in 3:k)
  {
    J=seq(2,(i-1)) 
    infect=sum(J*B[1,])*(1-pi)/pi # Estimated number infected
    alpha=max(beta*(N-(i-1)-infect)/N,0) # i-1 individuals have been removed
    p=alpha/(alpha+gamma)  # Probability an event is a birth
    q=1-p        # Probability an event is a death
    phi=exp(-(alpha+gamma)*t[i-1])
    psi=(1-q)*(1-phi)
    g=pi*q*(1-pi)*(alpha+gamma)*phi/(q-(q-pi)*phi)^2
    r=pi*(q+ (1-q)*phi)/(q-(q-pi)*phi)
    # pi update
    pi=(q-(q-pi)*phi)
    h=(1-pi-psi)/((1-pi)*(1-psi)) # Uses the updated pi
    #    print(c(h,alpha,q,i,infect))
    SS=sum(J*r^(J-1)*B[1,])  # E[R_{i-1} r^{R_{i-1}-1}]
    logl=logl+log(g)+log(SS) # Update log likelihood
    B=B%*%Mmat(h,r,i-1)
    B=B/sum(B)    # Update B
  }  
  J=seq(2,k) 
  infect=sum(J*B[1,])*(1-pi)/pi # Estimated number infected
  alpha=max(beta*(N-(i-1)-infect)/N,0)
  p=alpha/(alpha+gamma)
  q=1-p
  phi=exp(-(alpha+gamma)*tau)
  r=pi*(q+(1-q)*phi)/(q-(q-pi)*phi)
  SS=sum(r^J*B[1,])
  logl=logl+log(SS)
  logl
}
