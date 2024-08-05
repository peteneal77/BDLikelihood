#
# Exact GSE likelihood
# - Constructed and run Abikiliki data
# 
# Returns log-likelihood

library(expm)

#
# Qkk
#

MatQ1=function(alpha,gamma,N,k)
{
  Q=matrix(0,ncol=(N+1-k),nrow=(N+1-k))
  for(j in 1:(N+1-k)) Q[j,j]=-(alpha*(N-k+1-j)+N*gamma)*j/N
  for(j in 2:(N+1-k)) Q[(j-1),j]=alpha*(N-k+1-(j-1))*(j-1)/N
  Q
}


#
# Qk,k+1
#

MatQ2=function(gamma,N,k)
{
  Q=matrix(0,nrow=(N-k+1),ncol=(N-k))
  for(i in 2:(N-k+1)) Q[i,(i-1)]=i*gamma
  Q
}


Elike_GSEI=function(beta,gamma,t,st,N)
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
  
  # Above sets up inference based on removal times up to time st.
  
  k=length(t)+1
  u=rep(0,N)
  u[1]=1
  Q00=MatQ1(beta,gamma,N-1,0)
  Q01=MatQ2(gamma,N-1,0)
  LL=-u%*%solve(Q00)%*%Q01
  for(i in 1:(k-2))
  {
    D1=MatQ1(beta,gamma,N-1,i)
    D2=MatQ2(gamma,N-1,i)
    LL=LL%*%expm(D1*t[i])%*%D2
  }
  D1=MatQ1(beta,gamma,N-1,k-1)
  LL=LL%*%expm(D1*t[k-1])
  tilD2=diag(gamma*seq(1,N+1-k))
  D1=MatQ1(beta,gamma,N-1,k)
  nD1=length(D1[1,])+1
  tilD1=matrix(0,ncol=nD1,nrow=nD1)
  tilD1[2:nD1,2:nD1]=D1
  LL=LL%*%tilD2%*%expm(tau*tilD1)
  log(sum(LL))
}


