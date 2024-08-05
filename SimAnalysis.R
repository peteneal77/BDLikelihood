#
# Covid - Control effects estimation (Simulation study)
#
# Simulated data set: DailyDataS2.csv (first column day counter)
# Control Measures: SimControl2.csv

library(MASS) # Needed for simulating from multivariate normal

source("approxBD_functions.R") # Likelihood function
source("parameter_calc.R") # Parameter calculation function

source("daily_convert.R") # Sets and updates removal times
source("Covid_MCMC2.R") # MCMC algorithm

NC=5 # Number of countries
NT=6 # Number of control measures

R0_Country=c(4,3.5,4.5,3.8,4.2) # Initial R0 in each country
# Parameters for "DailyDataS2.csv": R0_Country=c(4,3.5,4.5,3.8,4.2)

theta_Control=c(0.9,0.8,0.8,0.75,0.35) # Effect of control measures
# Parameters for "DailyDataS2.csv":  theta_Control=c(0.9,0.8,0.8,0.75,0.35)

Country_Control=read.csv("SimControl2.csv",row.names = 1) 
Country_Control

DayZ=c(55,50,45,48,50) # Days of data after the first case observed
# Parameters for "DailyDataS2.csv":  DayZ=c(55,50,45,48,50) 

mu=0.1 # 1/mu - mean infectious period
d=0.1 # Proportion of cases detected 
# Parameters for "DailyDataS2.csv": mu=0.1, d=c(0.1,0.08,0.09,0.11,0.12)

beta_Country=R0_Country*mu
beta_Final=rep(0,NC) 
# beta_Final=c(-0.15, -0.01,  0.07, -0.24,0.17)

Tdaily1=read.csv("DailyDataS2.csv",row.names = 1) # Read in data
Tdaily=as.matrix(Tdaily1)
Nday=length(Tdaily1[,1]) # Total number of days in analysis

#
# Setup control measures by day. 1 - Control measure in place, 0 - control measure not in place.
#

Day_Control=array(1,dim=c(NC,NT,Nday)) 
# dimension - country, control measures, day
for(i in 1:NC)
{
  for(j in 1:NT)
  {
    Day_Control[i,j,1:min(Country_Control[i,j],Nday)]=0
  }
}

NT=NT-1 #subtract 1 for country effect.

#
# Setting parameters 
#

Nrun=10000
betaI=rep(0.4,5) # Initial betaI values corresponding to R0=4
contI=rep(0.7,5) # Initial NPI values (shared)
theta=c(betaI,mu,d,contI,beta_Final)
theta_fix=c(rep(1,5),0,0,rep(1,5),rep(1,5)) # Fix mu and d (same)
#theta_fix=rep(1,17) # No parameters fixed
theta_prob=c(rep(0,6),rep(1,6),rep(2,5)) 
# Identifies parameters; 0 - positive, 1 - probability, 2 - real
theta_int=rep(0,17) # No integer parameters
sigma=0.00001 # Initial value for Covariance matrix Sigma = sigma I

Nup=round(0.1*colSums(Tdaily1),0) # Update 10% of removal (death) times at each iteration.

# Set up priors
# mu ~ Gamma (10,100) ; d ~ Beta (10,90)
# NPIs ~ U(0,1)
# beta initial ~ Gamma (10,25); beta final ~ N(0,0.1^2)

prior=matrix(1,ncol=2,nrow=17)
prior[6,]=c(10,100)
prior[7,]=c(10,90)

for(i in 1:5) prior[i,]=c(10,25)
for(i in 13:17) prior[i,]=c(0,0.1)

# Run MCMC with initial values

A=Covid_MCMC(Nrun,Tdaily,DayZ,Country_Control,
                      theta,theta_fix,theta_prob,theta_int
                      ,prior,sigma,Nup)

# vary - mu and d varying, fix - mu and d fixed

write.csv(A,file="MCMC Output/EuroSimS2_fix.csv",row.names = FALSE)

theta=A[Nrun,]

A2=Covid_MCMC(Nrun,Tdaily,DayZ,Country_Control,
             theta,theta_fix,theta_prob,theta_int
             ,prior,sigma,Nup)

# vary - mu and d varying, fix - mu and d fixed

write.csv(A2,file="MCMC Output/EuroSimS2_fix2.csv",row.names = FALSE)

theta=A2[Nrun,]

Nrun=c(10000,10000,10000,50000)

A3=Covid_MCMC(Nrun,Tdaily,DayZ,Country_Control,
              theta,theta_fix,theta_prob,theta_int
              ,prior,sigma,Nup)

# vary - mu and d varying, fix - mu and d fixed

write.csv(A3,file="MCMC Output/EuroSimS2_fix3.csv",row.names = FALSE)


