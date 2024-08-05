#
# European Covid data analysis file
#
#
#
#
#
#

library(MASS) # Needed for simulating from multivariate normal

source("approxBD_functions.R") # Likelihood function
source("parameter_calc.R") # Parameter calculation function for birth, death and detection parameters

source("daily_convert.R") # Sets and updates removal times
source("Covid_MCMC.R") # MCMC algorithm

# Sets up requirements for running preparation file.

muD=18  # Mean R -> D
alphaD=0 # Shape parameter for R -> D Gamma (alphaD, alphaD/muD). alphaD = 0 fixed infectious period.
ChooCoun=seq(1:11)
ChCont=c(1,2,3,4,5,6,7)
# ChCont=c(2,5,6,7)

source("EuroCovidPrep.R") 
# Run preparation of data with chosen R->D distribution, countries and control measures.

NC=length(ChooCoun) # Number of countries chosen
NT=length(ChCont) # Number of control measures including final


R0_Country=rep(3,NC) # Initial R0 in each country
theta_Control=rep(0.8,(NT-1)) # Effect of control measures  minus individual country effect

beta_Final=rep(0,NC) # Effect of final intervention (relative) - set up for exponential - initially equal to 0

Country_Control=Control_Int # Control Measures

mu=2/13 # 1/mu - mean infectious period  - mu=0.1 corresponds to a mean of 10 days
d=0.1 # Proportion of cases detected [keep constant throughout]

beta_Country=R0_Country*mu # Set initial beta values

# Finalise daily counts data for analysis in the MCMC.

Tdaily1=EDeaths
Tdaily=as.matrix(Tdaily1)
Nday=length(Tdaily[,1])

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

NT=NT-1 # Minus 1 for the individual control measures

Nrun=10000 # Single run of 10000 iterations to initialise.
Tdata=Tdaily

theta=c(beta_Country,mu,d,theta_Control,beta_Final) # Initial parameter values

theta_fix=c(rep(1,NC),0,0,rep(1,NT),rep(1,NC)) # Parameters mu and d fix
#theta_fix=c(rep(1,NC),1,0,rep(1,NT),rep(1,NC)) # Fix d but not mu
theta_prob=c(rep(0,(NC+1)),rep(1,(NT+1)),rep(2,NC)) 
# Set parameter class 0 - positive, 1 - probability, 2 - real
theta_int=rep(0,length(theta)) # Sets integer parameters, there are none.
sigma=0.000001


EDeaths=as.matrix(EDeaths)
Nup=round(0.1*colSums(EDeaths),0)
# Update 10% of death (removal) times at each iteration.

# Set up priors
# Gamma (10,100) prior on mu and Beta (10,90) prior on d if these are updated.
# U(0,1) priors on NPIs
# Gamma (10,25) priors on beta initial and N(0,0.1^2 prior) on beta final
prior=matrix(1,ncol=2,nrow=length(theta))
prior[(NC+1),]=c(10,100)
prior[(NC+2),]=c(10,90)
prior[1:NC,1]=10
prior[1:NC,2]=25
for(ii in (length(theta)+1-NC):length(theta)) prior[ii,]=c(0,0.1) # N(0,0.1^2) prior on end effect

B=Covid_MCMC(Nrun,Tdaily,DayE,Country_Control,
              theta,theta_fix,theta_prob,theta_int
              ,prior,sigma,Nup)
# Countries (Euro full) F - Fixed, V - vary (mu first then d), k counter, abcdefg - NPIs
write.csv(B,file="Euro_FF1_65_1234567.csv",row.names = FALSE)


theta=0
for(i in 1:length(B[Nrun,])) theta[i]=B[Nrun,i]

B2=Covid_MCMC(Nrun,Tdaily,DayE,Country_Control,
             theta,theta_fix,theta_prob,theta_int
             ,prior,sigma,Nup)

write.csv(B2,file="Euro_FF2_65_1234567.csv",row.names = FALSE)

theta=B2[Nrun,]

Nrun=c(10000,10000,10000,50000)

B3=Covid_MCMC(Nrun,Tdaily,DayE,Country_Control,
              theta,theta_fix,theta_prob,theta_int
              ,prior,sigma,Nup)

write.csv(B3,file="Euro_FF3_65_1234567.csv",row.names = FALSE)
