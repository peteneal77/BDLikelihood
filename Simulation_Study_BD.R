#
# Comparison of birth-death likelihoods (exact and approximate)
#
#
#

source("BD_Exact_Like.R") # Exact log-likelihood
source("BD_Approx_Like.R") # Approximate log-likelihood
source("BD_Hybrid_Like.R") # Hybrid log-likelihood

source("BD_sim.R")

set.seed(73629231) # Allows repeatibility

N=200 # All simulations until 200 detected deaths.
LX=50 # All hybrid likelihoods use the exact likelihood for the first 50 detected deaths

#
# Homogeneous mixing with alpha=0.2, mu=0.1 and d=1
#

alpha1=rep(0.2,N)
mu1=rep(0.1,N)
d1=rep(1,N)

# Simulates data until a realisation with N=200 detected deaths occurs.

inter1=0
while(length(inter1)<(N-1))
{
  inter1=BD(N,alpha1,mu1,d1)
}

# write.csv(inter1,file="Sim_data1.csv",row.names = FALSE) # Data analysed

beta1=seq(0.1,0.3,0.001)
gamma1=seq(0.01,0.20,0.001)  # Grid points for the likelihood to be evaluated at.

E1=matrix(0,nrow=length(beta1),ncol=length(gamma1))
A1=matrix(0,nrow=length(beta1),ncol=length(gamma1))
H1=matrix(0,nrow=length(beta1),ncol=length(gamma1))

a1_before=proc.time()

for(i in 1:length(beta1))
{
  for(j in 1:length(gamma1))
  {
    A1[i,j]=like_BD_approx(rep(beta1[i],N),rep(gamma1[j],N),d1,inter1)
  }
}

a1_after=proc.time()

a1_after-a1_before

#user  system elapsed 
#39.51    0.03   39.63 

filled.contour(beta1,gamma1,exp(A1),xlab=expression(alpha),ylab=expression(mu),main="Approximate")

# Contour plot on grayscale

gg=gray.colors(n=23,start=1,end=0,gamma=1)

filled.contour(beta1,gamma1,exp(A1),xlab=expression(alpha),ylab=expression(mu),
               levels=seq(0,2.2*10^78,10^77),col=gg,main="Approximate")

#
#
#

e1_before=proc.time()

for(i in 1:length(beta1))
{
  print(i)
  for(j in 1:length(gamma1))
  {
    E1[i,j]=like_BD_exact(rep(beta1[i],N),rep(gamma1[j],N),d1,inter1)
  }
}

e1_after=proc.time()

e1_after-e1_before

filled.contour(beta1,gamma1,exp(E1),xlab=expression(alpha),ylab=expression(mu),main="Exact") 

#> e1_after-e1_before
# user   system  elapsed 
# 22985.06    19.71 23064.00 

# Pre-calculated exact likelihoods
#E1=read.csv("Output/Sim_data1_homo_exact.csv")
#E1=as.matrix(E1)

gg=gray.colors(n=23,start=1,end=0,gamma=1)

filled.contour(beta1,gamma1,exp(E1),xlab=expression(alpha),ylab=expression(mu),
               levels=seq(0,2.2*10^78,10^77),col=gg,main="Exact")

h1_before=proc.time()

for(i in 1:length(beta1))
{
  print(i)
  for(j in 1:length(gamma1))
  {
    H1[i,j]=like_BD_hybrid(rep(beta1[i],N),rep(gamma1[j],N),d1,inter1,LX)
  }
}

h1_after=proc.time()

h1_after-h1_before

# > h1_after-h1_before
#    user  system elapsed 
# 493.24    1.56  496.92 

filled.contour(beta1,gamma1,exp(H1),xlab=expression(alpha),ylab=expression(mu),main="Hybrid") 

# Pre-calculated exact likelihoods
H1=read.csv("Output/Sim_data1_homo_hybrid.csv")
H1=as.matrix(H1)

# Contour plot on grayscale

gg=gray.colors(n=23,start=1,end=0,gamma=1)

filled.contour(beta1,gamma1,exp(H1),xlab=expression(alpha),ylab=expression(mu),
               levels=seq(0,2.2*10^78,10^77),col=gg,main="Hybrid")


D1=E1
for(i in 1:length(beta1))
{
  print(i)
  for(j in 1:length(gamma1))
  {
    D1[i,j]=min(E1[i,j],A1[i,j])-max(E1[i,j],A1[i,j]) 
  }  
}

filled.contour(beta1,gamma1,exp(D1),xlab=expression(alpha),ylab=expression(mu),main="Difference") 

# Contour plot on grayscale

gg=gray.colors(n=26,start=1,end=0,gamma=1)

filled.contour(beta1,gamma1,exp(D1),xlab=expression(alpha),ylab=expression(mu),
               levels=seq(0,1,1/25),col=gg,main="Difference")


#
# 01.03.24 - Simulation study - calculates the hybrid likelihood at each possible value of M.
#

beta1=seq(0.1,0.3,0.001)
gamma1=seq(0.01,0.20,0.001) 

Hfull=array(0,dim=c(length(beta1),length(gamma1),(N-1)))

hfull_before=proc.time()

for(i in 1:length(beta1))
{
  print(i)
  for(j in 1:length(gamma1))
  {
    ll=like_BD_hybrid_ext(rep(beta1[i],N),rep(gamma1[j],N),d1,inter1)
    Hfull[i,j,]=ll[2:N]
  }
}

hfull_after=proc.time()

hfull_after-hfull_before

# write.csv(Hfull,file="Output/Sim_data1_homo_combo.csv",row.names = FALSE)

# XX=read.csv("Output/Sim_data1_homo_combo.csv")

#for(SEED in 2:N)
#{  
#  XA=Hfull[,,(SEED-1)]
#  write.csv(XA,file=paste0("Output/Sim_data1_homo_combo_",SEED,".csv"),row.names = FALSE)
#  print(SEED)
#}


# Reads in the data
# 
#

N=200

AA=array(0,dim=c(201,191,199))
for(SEED in 2:N)
{  
  XB=read.csv(file=paste0("Output/Sim_data1_homo_combo_",SEED,".csv"))
  XB=as.matrix(XB)
  AA[,,(SEED-1)]=XB
  print(SEED)
}

Hfull=AA

tt=131
plot.ts(exp(Hfull[,tt,199]-max(Hfull[,tt,])),type="l")
for(i in 2:198) lines(exp(Hfull[,tt,i]-max(Hfull[,tt,])))

ta=131
plot.ts(exp(Hfull[ta,,199]-max(Hfull[ta,,])),type="l")
for(i in 2:198) lines(exp(Hfull[ta,,i]-max(Hfull[ta,,])))

tt=131
plot.ts(exp(Hfull[ta,,199]),type="l")
for(i in 2:198) lines(exp(Hfull[ta,,i]))

plot.ts(Hfull[,tt,199])
for(i in 2:198) lines(Hfull[,tt,i],col=i)
