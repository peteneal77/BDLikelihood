#
# Simulated difference
#


source("BD_Exact_Like.R") # Exact log-likelihood
source("BD_Approx_Like.R") # Approximate log-likelihood
source("BD_Hybrid_Like.R") # Hybrid log-likelihood

source("BD_sim.R")

set.seed(938777113)

MN=100000

N=5

alpha1=rep(0.4,N)
mu1=rep(0.1,N)
d1=rep(1,N)

EL=0
AL=0

for(i in 1:MN)
{
  inter1=0
  while(length(inter1)<(N-1))
  {
    inter1=BD(N,alpha1,mu1,d1)
  }
  EL[i]=like_BD_exact(alpha1,mu1,d1,inter1)
  AL[i]=like_BD_approx(alpha1,mu1,d1,inter1)
}
