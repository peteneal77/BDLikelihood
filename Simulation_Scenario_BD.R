#
# Comparison of birth-death likelihoods (exact and approximate)
#
#
#

source("BD_Sim.R") # Simulation code

source("MCMC_BD.R")  # MCMC code 
# all sub-functions likelihoods, parameters and priors loaded in from MCMC_BD.R


N=200 # All simulations until 200 detected deaths.
# LX=50 # All hybrid likelihoods use the exact likelihood for the first 50 detected deaths
# Note that like_BD_hybrid_Fix fixes LX (L)=50

Nrun=c(5000,5000,5000,20000) # MCMC runs - 3 * 5000 burn-ins to fix proposal variance
sigtuna=0.01 
# Initial standard deviation for each component of random walk except any discrete variables which have set initial standard deviation 1.

# Scenario 0
#
# Time-homogeneous with alpha=0.2, mu=0.1 and d=0.1
#
# Simulation_Study_BD.R
#


# Scenario 1
#
# Time-homogeneous with changepoint after detected death LC=50: mu=0.1 and d=0.4
# Before: alpha=0.3; After: alpha=0.05

set.seed(62527592)

theta2=c(0.3,0.05,0.1,0.4,50) # (alpha1, alpha2, mu, d, LC)

Pa2=para2(theta2,N) # Generate alpha, mu, d vectors
alpha2=Pa2$alpha
mu2=Pa2$mu
d2=Pa2$d

# Simulate inter-arrival times.
# Repeat simulation until we get N=200 detected deaths

inter2=0
while(length(inter2)<(N-1))
{
  inter2=BD(N,alpha2,mu2,d2)
}

# 

theta_type2=c(1,1,1,2,3) # Parameter types 1-positive, 2-probability, 3-discrete (integer)
theta_fix2=c(1,1,1,1,1) # Parameters to update 0 -fixed, 1-update
theta_prior2=matrix(0,ncol=3,nrow=5) # Set priors
theta_prior2[1,]=c(2,5,5/theta2[1])
theta_prior2[2,]=c(2,5,5/theta2[2])
theta_prior2[3,]=c(2,5,5/theta2[3])
theta_prior2[4,]=c(3,10*theta2[4],10*(1-theta2[4]))
theta_prior2[5,]=c(1,1,200)

theta_initial2=c(0.25,0.1,0.2,0.7,100) # INitial values for MCMC run

# Run MCMC

QE2=MCMC_BD(inter2,like_BD_exact,para2,theta_initial2,Nrun,sigtuna,theta_type2,theta_prior2,theta_fix2)
QA2=MCMC_BD(inter2,like_BD_approx,para2,theta_initial2,Nrun,sigtuna,theta_type2,theta_prior2,theta_fix2)
QH2=MCMC_BD(inter2,like_BD_hybrid_Fix,para2,theta_initial2,Nrun,sigtuna,theta_type2,theta_prior2,theta_fix2)

#write.csv(QE2,file="QE2.csv",row.names = FALSE)
#write.csv(QA2,file="QA2.csv",row.names = FALSE)
#write.csv(QH2,file="QH2.csv",row.names = FALSE)


# Scenario 2
#
# Time-inhomogeneous birth rate with decay

set.seed(54341028)

theta3=c(0.25,0.2,0.1,0.25) # (beta,phi,gamma,delta)

Pa3=para3(theta3,N) # Generate alpha, mu, d vectors
alpha3=Pa3$alpha
mu3=Pa3$mu
d3=Pa3$d

# Simulate inter-arrival times.
# Repeat simulation until we get N=200 detected deaths

inter3=0
while(length(inter3)<(N-1))
{
  inter3=BD(N,alpha3,mu3,d3)
}

# 

theta_type3=c(1,1,1,2) # Parameter types 1-positive, 2-probability, 3-discrete (integer)
theta_fix3=c(1,1,1,1) # Parameters to update 0 -fixed, 1-update
theta_prior3=matrix(0,ncol=3,nrow=4) # Set priors
theta_prior3[1,]=c(2,5,5/theta3[1])
theta_prior3[2,]=c(3,5*theta3[2],5*(1-theta3[2]))
theta_prior3[3,]=c(2,5,5/theta3[3])
theta_prior3[4,]=c(3,5*theta3[4],5*(1-theta3[4]))

theta_initial3=c(0.2,0.5,0.2,0.6)

# Run MCMC

QE3=MCMC_BD(inter3,like_BD_exact,para3,theta_initial3,Nrun,sigtuna,theta_type3,theta_prior3,theta_fix3)
QA3=MCMC_BD(inter3,like_BD_approx,para3,theta_initial3,Nrun,sigtuna,theta_type3,theta_prior3,theta_fix3)
QH3=MCMC_BD(inter3,like_BD_hybrid_Fix,para3,theta_initial3,Nrun,sigtuna,theta_type3,theta_prior3,theta_fix3)

#write.csv(QE3,file="QE3.csv",row.names = FALSE)
#write.csv(QA3,file="QA3.csv",row.names = FALSE)
#write.csv(QH3,file="QH3.csv",row.names = FALSE)

# Scenario 3
#
# Time-homogeneous model with changepoint

set.seed(30720956)

theta4=c(0.2,0.1,0.2,0.15,0.15,0.5,100) # (beta1,gamma1,delta1,beta2,gamma2,delta2)

Pa4=para4(theta4,N) # Generate alpha, mu, d vectors
alpha4=Pa4$alpha
mu4=Pa4$mu
d4=Pa4$d

# Simulate inter-arrival times.
# Repeat simulation until we get N=200 detected deaths

inter4=0
while(length(inter4)<(N-1))
{
  inter4=BD(N,alpha4,mu4,d4)
}

# 

theta_type4=c(1,1,2,1,1,2,3) 
# Parameter types 1-positive, 2-probability, 3-discrete (integer)
theta_fix4=rep(1,7) # Parameters to update 0 -fixed, 1-update
theta_prior4=matrix(0,ncol=3,nrow=7) # Set priors
theta_prior4[1,]=c(2,5,5/theta4[1])
theta_prior4[2,]=c(2,5,5/theta4[2])
theta_prior4[3,]=c(3,5*theta4[3],5*(1-theta4[3]))
theta_prior4[4,]=c(2,5,5/theta4[4])
theta_prior4[5,]=c(2,5,5/theta4[5])
theta_prior4[6,]=c(3,5*theta4[6],5*(1-theta4[6]))
theta_prior4[7,]=c(1,1,N)

theta_initial4=c(0.15,0.05,0.4,0.1,0.2,0.2,50)

# Run MCMC

QE4=MCMC_BD(inter4,like_BD_exact,para4,theta_initial4,Nrun,sigtuna,theta_type4,theta_prior4,theta_fix4)
QA4=MCMC_BD(inter4,like_BD_approx,para4,theta_initial4,Nrun,sigtuna,theta_type4,theta_prior4,theta_fix4)
QH4=MCMC_BD(inter4,like_BD_hybrid_Fix,para4,theta_initial4,Nrun,sigtuna,theta_type4,theta_prior4,theta_fix4)

#write.csv(QE4,file="QE4.csv",row.names = FALSE)
#write.csv(QA4,file="QA4.csv",row.names = FALSE)
#write.csv(QH4,file="QH4.csv",row.names = FALSE)

# Scenario 4
#
# Time-inhomogeneous model with changing rates in all parameters
# beta decreasing, gamma and delta increasing

set.seed(46208046)

theta5=c(0.4,0.25,0.1,0.5,0.8,0.25) 
# (beta,phi_beta,gamma,phi_gamma,delta,phi_delta)

Pa5=para5(theta5,N) # Generate alpha, mu, d vectors
alpha5=Pa5$alpha
mu5=Pa5$mu
d5=Pa5$d

# Simulate inter-arrival times.
# Repeat simulation until we get N=200 detected deaths

inter5=0
while(length(inter5)<(N-1))
{
  inter5=BD(N,alpha5,mu5,d5)
}

# 

theta_type5=c(1,2,1,2,1,2) 
# Parameter types 1-positive, 2-probability, 3-discrete (integer)
theta_fix5=rep(1,6) # Parameters to update 0 -fixed, 1-update
theta_prior5=matrix(0,ncol=3,nrow=6) # Set priors
theta_prior5[1,]=c(2,5,5/theta5[1])
theta_prior5[2,]=c(3,5*theta5[2],5*(1-theta5[2]))
theta_prior5[3,]=c(2,5,5/theta5[3])
theta_prior5[4,]=c(3,5*theta5[4],5*(1-theta5[4]))
theta_prior5[5,]=c(2,5,5/theta5[5])
theta_prior5[6,]=c(3,5*theta5[6],5*(1-theta5[6]))

theta_initial5=c(0.15,0.5,0.05,0.333,0.6,0.5)

# Run MCMC

QE5=MCMC_BD(inter5,like_BD_exact,para5,theta_initial5,Nrun,sigtuna,theta_type5,theta_prior5,theta_fix5)
QA5=MCMC_BD(inter5,like_BD_approx,para5,theta_initial5,Nrun,sigtuna,theta_type5,theta_prior5,theta_fix5)
QH5=MCMC_BD(inter5,like_BD_hybrid_Fix,para5,theta_initial5,Nrun,sigtuna,theta_type5,theta_prior5,theta_fix5)

#write.csv(QE5,file="QE5.csv",row.names = FALSE)
#write.csv(QA5,file="QA5.csv",row.names = FALSE)
#write.csv(QH5,file="QH5.csv",row.names = FALSE)

# Scenario 3 - Fixed L=100
#
# Time-homogeneous model with changepoint

set.seed(30720956)

theta4=c(0.2,0.1,0.2,0.15,0.15,0.5,100) # (beta1,gamma1,delta1,beta2,gamma2,delta2)

Pa4=para4(theta4,N) # Generate alpha, mu, d vectors
alpha4=Pa4$alpha
mu4=Pa4$mu
d4=Pa4$d

# Simulate inter-arrival times.
# Repeat simulation until we get N=200 detected deaths

inter4=0
while(length(inter4)<(N-1))
{
  inter4=BD(N,alpha4,mu4,d4)
}

# 

theta_type4=c(1,1,2,1,1,2,3) 
# Parameter types 1-positive, 2-probability, 3-discrete (integer)
theta_fix4=c(rep(1,6),0) # Parameters to update 0 -fixed, 1-update
theta_prior4=matrix(0,ncol=3,nrow=7) # Set priors
theta_prior4[1,]=c(2,5,5/theta4[1])
theta_prior4[2,]=c(2,5,5/theta4[2])
theta_prior4[3,]=c(3,5*theta4[3],5*(1-theta4[3]))
theta_prior4[4,]=c(2,5,5/theta4[4])
theta_prior4[5,]=c(2,5,5/theta4[5])
theta_prior4[6,]=c(3,5*theta4[6],5*(1-theta4[6]))
theta_prior4[7,]=c(1,1,N)

theta_initial4=c(0.15,0.05,0.4,0.1,0.2,0.2,100) # Fix L=100

# Run MCMC

QE4b=MCMC_BD(inter4,like_BD_exact,para4,theta_initial4,Nrun,sigtuna,theta_type4,theta_prior4,theta_fix4)
QA4b=MCMC_BD(inter4,like_BD_approx,para4,theta_initial4,Nrun,sigtuna,theta_type4,theta_prior4,theta_fix4)
QH4b=MCMC_BD(inter4,like_BD_hybrid_Fix,para4,theta_initial4,Nrun,sigtuna,theta_type4,theta_prior4,theta_fix4)

#write.csv(QE4b,file="QE4b.csv",row.names = FALSE)
#write.csv(QA4b,file="QA4b.csv",row.names = FALSE)
#write.csv(QH4b,file="QH4b.csv",row.names = FALSE)
