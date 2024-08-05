#
# Code for simulating Covid-19 type data
# - Multiple outbreaks in different countries
# - Shared NPI effects across countries
# - Shared infectious period distribution, Exp(mu)
# - Country R_0, death probability d and final control measure effect.
# - Data generated with no delay R->D (equivalent to a fixed length period)


#
#
#

NC=5 # Number of countries

R0_Country=c(4,3.5,4.5,3.8,4.2) # Initial R0 in each country
# Parameters for "DailyDataS2.csv": R0_Country=c(4,3.5,4.5,3.8,4.2)

theta_Control=c(0.9,0.8,0.8,0.75,0.35) # Effect of control measures
# Parameters for "DailyDataS2.csv":  theta_Control=c(0.9,0.8,0.8,0.75,0.35)

Country_Control=read.csv("SimControl2.csv",row.names = 1) # Control measures
Country_Control

#
# Control measures after a given day 
#

source("simulation_P_daily.R")

DayZ=c(55,50,45,48,50) # Days of data after the first case observed

mu=0.1 # 1/mu - mean infectious period
d=c(0.1,0.08,0.09,0.11,0.12) # Proportion of cases detected [keep constant throughout]
# Parameters for "DailyDataS2.csv": mu=0.1, d=0.1 

beta_Final=c(-0.15, -0.01,  0.07, -0.24,0.17) # Final effect in each country.

Tdaily=matrix(0,nrow=max(DayZ),ncol=NC) # A matrix to store daily counts from different countries

#
# theta_gen_daily calculates the daily infection rate
#

theta_gen_daily=function(Pday,R0,theta_Control,Day_Control,mu)
{
  betaday=rep(mu*R0,Pday) # Initial rate - parameterise as R0 and recovery rate
  DayA=matrix(1,ncol=length(Day_Control),nrow=Pday)
  for(j in 1:length(Day_Control)) DayA[1:min(Day_Control[j],Pday),j]=0
  # Allows for control measures introduced after period
  for(i in 1:Pday)
  {
    if(sum(DayA[i,])>0)
    {
      betaday[i]=betaday[i]*prod(theta_Control[DayA[i,]==1])
    }
  }
  betaday
}


for(i in 1:NC)
{
  print(i)
  gammaday=rep(mu,DayZ[i]) # No change in lifetime
  detpday=rep(d[i],DayZ[i])  # No change in detection
  betaday=theta_gen_daily(DayZ[i],R0_Country[i],c(theta_Control,exp(beta_Final[i]))
                          ,Country_Control[i,],mu)
  trix=0
  while(trix==0)
  {
    A=Sim_P_daily(DayZ[i],betaday,gammaday,detpday)
    trix=A$daily[DayZ[i]]
  }
  Tdaily[1:DayZ[i],i]=A$daily
}

colSums(Tdaily)

#write.csv(Tdaily,file="DailyDataS2.csv")

# 1740 1067 6785  345 4123 - Total number of cases




